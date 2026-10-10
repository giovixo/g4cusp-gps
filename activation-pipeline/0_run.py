"""
0_run.py
========
Step 0 of the activation pipeline: run the Geant4 activation simulation for a
set of monochromatic primary energies and produce one score CSV per energy.

For every energy the script
  1. creates an isolated work directory  <outdir>/work/<E>_MeV/  and links into
     it every file of the mass-model directory (GDML + included XML), so that
     runs never see each other's output and stale files cannot be picked up;
  2. writes the GPS macro (kept in the work directory for provenance);
  3. runs the Geant4 executable there, logging stdout/stderr to run.log;
  4. merges the per-thread ntuple files (<prefix>_run<N>_nt_<ntuple>_t<T>.csv)
     into  <outdir>/scorefile_<E>_MeV_<nprim>.csv  (plain CSV with a header row,
     sorted by EventID).  An energy with no recorded nuclide gives a header-only
     file, so "no activation" is distinguishable from "not simulated".
  5. records the run in  <outdir>/run_info.json  (source geometry, beam area,
     primaries, files, timing).  Later pipeline steps should read the source
     geometry from this file, so that the normalisation stays consistent.

Source geometry and normalisation
---------------------------------
Primaries start on a sphere of radius R (centre C) and point inwards with a
cosine-law angular distribution restricted to theta <= thetamax.  For an
omnidirectional integral flux F [p/cm^2/s] the number of primaries that
corresponds to an exposure time dt is

    N = dt * F * pi * R^2 * sin^2(thetamax)          (beam_area_cm2 = pi R^2 sin^2 thetamax)

This holds only if (a) the whole mass model is inside the sphere of radius
R*sin(thetamax) around C, and (b) the source sphere is inside the Geant4 world
volume.  Geant4 does NOT stop primaries that start outside the world: it moves
them straight across it and they miss the payload.  For the CUSP model (world
box 50 cm, payload bounding box x -49..47, y -70..40, z -46..142 mm) the
defaults use R = 13 cm centred on the payload, with the full cosine law
(thetamax = 90 deg).

Configuration
-------------
All settings have defaults for the CUSP model and can be set in a TOML file
(--config, see config_cusp.toml) and overridden on the command line.
Relative paths in the config file are resolved relative to the config file.

Usage
-----
    python 0_run.py                                  # CUSP defaults
    python 0_run.py --config config_cusp.toml
    python 0_run.py --energies 10 30 100 --nprim 1e6 --threads 8
    python 0_run.py --g4-args=-n --outdir output_nact   # neutron activation on
    python 0_run.py --dry-run                        # write macros, print commands

Completed energies are skipped on rerun (same nprim and output file present),
unless --overwrite is given.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Defaults (CUSP mass model)
# ---------------------------------------------------------------------------

DEFAULTS = {
    # Geant4 application
    "executable":    str(HERE.parent / "Debug" / "cusp-activation"),
    "geometry_dir":  str(HERE.parent / "gdml-mass-model"),
    "threads":       0,          # 0 = all cores (passed as -t)
    "g4_args":       [],         # extra command-line arguments, e.g. ["-n"]
    "pre_commands":  ["/cusp/stepping/verbose 0"],   # macro commands before /run/beamOn
    "thread_file_regex": r"_run(?P<run>\d+)_nt_(?P<ntuple>\w+?)_t(?P<thread>\d+)\.csv$",
    "ntuple":        "Events",
    # Primary source (GPS)
    "particle":      "proton",
    "energies":      [5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400, 700],
    "nprim":         100_000_000,
    "radius_cm":     13.0,
    "centre_cm":     [0.0, -1.5, 4.8],
    "thetamax_deg":  90.0,
    # Output
    "outdir":        "output",
    "keep_work":     True,       # keep the per-energy work directories (logs, macros, thread files)
}

PATH_KEYS = ("executable", "geometry_dir", "outdir")
FALLBACK_COLUMNS = ["EventID", "Isotope", "Lifetime", "Volume"]


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def load_config(path: Path) -> dict:
    """Read a TOML config file; resolve relative paths against its directory."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    unknown = set(cfg) - set(DEFAULTS)
    if unknown:
        sys.exit(f"ERROR: unknown key(s) in {path}: {sorted(unknown)}")
    for key in PATH_KEYS:
        if key in cfg and not Path(cfg[key]).expanduser().is_absolute():
            cfg[key] = str((path.parent / cfg[key]).resolve())
    return cfg


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run the Geant4 activation simulation for a set of monochromatic energies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Configuration", 1)[1],
    )
    p.add_argument("--config", type=Path, help="TOML configuration file")
    p.add_argument("--executable", help="Geant4 executable")
    p.add_argument("--geometry-dir", dest="geometry_dir",
                   help="directory with the GDML mass model (all files are linked into each work dir)")
    p.add_argument("--threads", "-t", type=int, help="worker threads, passed as -t (0 = all cores)")
    p.add_argument("--g4-args", dest="g4_args", type=shlex.split,
                   help='extra arguments for the executable, as one string, e.g. --g4-args="-n"')
    p.add_argument("--pre-command", dest="pre_commands", action="append",
                   help="macro command run before /run/beamOn (repeatable; replaces the configured list)")
    p.add_argument("--particle", help="GPS particle name")
    p.add_argument("--energies", nargs="+", type=float, help="primary energies [MeV]")
    p.add_argument("--nprim", type=float, help="primaries per energy (e.g. 1e8)")
    p.add_argument("--radius", dest="radius_cm", type=float, help="source sphere radius [cm]")
    p.add_argument("--centre", dest="centre_cm", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="source sphere centre [cm]")
    p.add_argument("--thetamax", dest="thetamax_deg", type=float,
                   help="maximum angle of the cosine-law distribution [deg]")
    p.add_argument("--outdir", help="output directory")
    p.add_argument("--no-keep-work", dest="keep_work", action="store_false", default=None,
                   help="delete the per-energy work directories after merging")
    p.add_argument("--overwrite", action="store_true", help="rerun energies that are already done")
    p.add_argument("--dry-run", action="store_true", help="write the macros and print the commands only")
    return p


def resolve_settings(args: argparse.Namespace) -> dict:
    cfg = dict(DEFAULTS)
    if args.config:
        cfg.update(load_config(args.config.resolve()))
    for key in DEFAULTS:
        val = getattr(args, key, None)
        if val is not None:
            cfg[key] = val
    cfg["nprim"] = int(cfg["nprim"])
    cfg["energies"] = sorted(float(e) for e in cfg["energies"])
    for key in PATH_KEYS:
        cfg[key] = str(Path(cfg[key]).expanduser().resolve())
    return cfg


def beam_area_cm2(radius_cm: float, thetamax_deg: float) -> float:
    return math.pi * radius_cm**2 * math.sin(math.radians(thetamax_deg))**2


# ---------------------------------------------------------------------------
# Per-energy run
# ---------------------------------------------------------------------------

def energy_tag(energy: float) -> str:
    """10.0 -> '10', 4.5 -> '4.5' (matches the parser's '<E>_MeV' pattern)."""
    return f"{energy:g}"


def score_filename(energy: float, nprim: int) -> str:
    return f"scorefile_{energy_tag(energy)}_MeV_{nprim:.0e}.csv"


def write_macro(path: Path, cfg: dict, energy: float) -> None:
    cx, cy, cz = cfg["centre_cm"]
    lines = list(cfg["pre_commands"]) + [
        "/gps/pos/type Surface",
        "/gps/pos/shape Sphere",
        f"/gps/pos/centre {cx:g} {cy:g} {cz:g} cm",
        f"/gps/pos/radius {cfg['radius_cm']:g} cm",
        "/gps/ang/type cos",
        "/gps/ang/mintheta 0 deg",
        f"/gps/ang/maxtheta {cfg['thetamax_deg']:g} deg",
        f"/gps/particle {cfg['particle']}",
        f"/gps/energy {energy:g} MeV",
        f"/run/beamOn {cfg['nprim']:d}",
    ]
    path.write_text("\n".join(lines) + "\n")


def prepare_workdir(workdir: Path, geometry_dir: Path) -> None:
    """Fresh work dir with symlinks to every file of the mass model."""
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    for f in geometry_dir.iterdir():
        if f.is_file():
            (workdir / f.name).symlink_to(f)


def read_thread_file(path: Path) -> tuple[list[str], list[list[str]]]:
    """Return (column names, data rows) of a Geant4 CSV ntuple file."""
    columns, rows = [], []
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("#"):
                parts = line.split()
                if len(parts) == 3 and parts[0] == "#column":
                    columns.append(parts[2])
            elif line.strip():
                rows.append(next(csv.reader([line])))
    return columns, rows


def merge_thread_files(workdir: Path, cfg: dict, dest: Path) -> dict:
    """Merge the thread files of one run into dest. Returns summary info."""
    pattern = re.compile(cfg["thread_file_regex"])
    files = []
    for f in sorted(workdir.iterdir()):
        m = pattern.search(f.name)
        if m and m.group("ntuple") == cfg["ntuple"]:
            files.append((int(m.group("run")), f))
    runs = {r for r, _ in files}
    if len(runs) > 1:
        raise RuntimeError(f"thread files from several runs in {workdir}: {sorted(runs)}")

    columns, rows = None, []
    for _, f in files:
        cols, data = read_thread_file(f)
        if columns is None:
            columns = cols
        elif cols != columns:
            raise RuntimeError(f"column mismatch in {f.name}: {cols} != {columns}")
        rows.extend(data)
    columns = columns or FALLBACK_COLUMNS

    if "EventID" in columns:
        i = columns.index("EventID")
        rows.sort(key=lambda r: int(float(r[i])))
        n_events = len({r[i] for r in rows})
    else:
        n_events = None

    tmp = dest.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(rows)
    tmp.replace(dest)
    return {"n_thread_files": len(files), "n_records": len(rows), "n_events_with_records": n_events}


def run_energy(cfg: dict, energy: float, outdir: Path, dry_run: bool) -> dict:
    tag = energy_tag(energy)
    workdir = outdir / "work" / f"{tag}_MeV"
    prepare_workdir(workdir, Path(cfg["geometry_dir"]))
    macro = workdir / "run.mac"
    write_macro(macro, cfg, energy)

    cmd = [cfg["executable"], "-t", str(cfg["threads"]), *cfg["g4_args"], macro.name]
    print(f"  cwd: {workdir}\n  cmd: {' '.join(cmd)}")
    if dry_run:
        return {}

    t0 = time.monotonic()
    with open(workdir / "run.log", "w") as log:
        proc = subprocess.run(cmd, cwd=workdir, stdout=log, stderr=subprocess.STDOUT)
    wall = time.monotonic() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Geant4 exited with code {proc.returncode} (see {workdir / 'run.log'})")

    dest = outdir / score_filename(energy, cfg["nprim"])
    info = merge_thread_files(workdir, cfg, dest)
    if not cfg["keep_work"]:
        shutil.rmtree(workdir)
    info.update({"energy_MeV": energy, "nprim": cfg["nprim"], "file": dest.name,
                 "wall_time_s": round(wall, 1),
                 "finished": dt.datetime.now().isoformat(timespec="seconds")})
    return info


# ---------------------------------------------------------------------------
# Run bookkeeping
# ---------------------------------------------------------------------------

def load_manifest(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_manifest(path: Path, manifest: dict) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    tmp.replace(path)


def source_block(cfg: dict) -> dict:
    return {
        "particle":      cfg["particle"],
        "shape":         "sphere surface, inward cosine law",
        "radius_cm":     cfg["radius_cm"],
        "centre_cm":     list(cfg["centre_cm"]),
        "thetamax_deg":  cfg["thetamax_deg"],
        "beam_area_cm2": beam_area_cm2(cfg["radius_cm"], cfg["thetamax_deg"]),
    }


def main(argv: list[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    cfg = resolve_settings(args)

    exe, geo = Path(cfg["executable"]), Path(cfg["geometry_dir"])
    if not args.dry_run and not os.access(exe, os.X_OK):
        sys.exit(f"ERROR: executable not found or not executable: {exe}")
    if not geo.is_dir() or not any(geo.glob("*.gdml")):
        sys.exit(f"ERROR: no GDML file in geometry directory {geo}")

    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "run_info.json"
    manifest = load_manifest(manifest_path)

    source = source_block(cfg)
    if manifest.get("source") and manifest["source"] != source:
        sys.exit(f"ERROR: source geometry differs from the one recorded in {manifest_path}.\n"
                 f"  recorded: {manifest['source']}\n  current:  {source}\n"
                 f"Use another --outdir to keep the normalisation of each output directory consistent.")
    manifest["source"] = source
    manifest.setdefault("runs", {})
    manifest["geant4"] = {"executable": str(exe), "threads": cfg["threads"], "g4_args": cfg["g4_args"],
                          "pre_commands": cfg["pre_commands"], "geometry_dir": str(geo)}

    print(f"Source: R = {cfg['radius_cm']:g} cm, centre = {cfg['centre_cm']} cm, "
          f"thetamax = {cfg['thetamax_deg']:g} deg, beam area = {source['beam_area_cm2']:.4g} cm^2")
    print(f"Energies [MeV]: {[energy_tag(e) for e in cfg['energies']]}, {cfg['nprim']:.3g} primaries each")

    failed = []
    for energy in cfg["energies"]:
        tag = energy_tag(energy)
        print(f"\n{'=' * 60}\n{cfg['particle']} {tag} MeV, {cfg['nprim']:.3g} primaries\n{'=' * 60}")

        done = manifest["runs"].get(tag)
        if (done and not args.overwrite and done.get("nprim") == cfg["nprim"]
                and (outdir / done["file"]).exists()):
            print(f"  already done ({done['file']}), skipped (use --overwrite to rerun)")
            continue

        try:
            info = run_energy(cfg, energy, outdir, args.dry_run)
        except (RuntimeError, OSError) as exc:
            print(f"  ERROR: {exc}")
            failed.append(tag)
            continue
        if args.dry_run:
            continue

        manifest["runs"][tag] = info
        save_manifest(manifest_path, manifest)
        print(f"  {info['n_records']} nuclides in {info['n_events_with_records']} events "
              f"({info['wall_time_s']} s) -> {info['file']}")

    if not args.dry_run:
        save_manifest(manifest_path, manifest)
        print(f"\nRun summary written to {manifest_path}")
    if failed:
        print(f"FAILED energies [MeV]: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
