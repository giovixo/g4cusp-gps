"""
run_postactivation.py
=====================
Step 4 of the activation pipeline: run the Geant4 post-activation program
(cusp-postactivation), which decays each active nuclide in each volume and
records the energy deposited in the scintillators, for every (volume, isotope)
pair of active_isotopes.pkl.

Allocation of the decays
------------------------
The background is the sum over pairs p of  A_p * f_p,  with A_p the
steady-state out-of-belt mean activity [Bq] and f_p the detector response per
decay.  With N_p simulated decays, the variance of the sum is
sum_p A_p^2 s_p^2 / N_p (s_p^2 = variance of one decay), which for a fixed
total sum N_p is minimised by N_p proportional to A_p s_p.  The weight used
here is the activity alone (s_p is unknown before the run), clipped to
[nmin, nmax] decays per pair and scaled to a total of about `budget`.

Weight of a pair: the same lag-weighted average as average_spectrum.py,

    A_p = sum_k w_k K_k[a_p],   a_p(tau) = sum_j N_j a_{j,p}(tau)

with a_{j,p} the activity per primary (activities.pkl), N_j the primaries of
energy E_j in one SPENVIS time step (accumulate_spectra.primaries_per_step),
K_k the step-averaged kernel and w_k the lag weights of the out-of-belt steps
(activation_history.lag_weights).  The sum of the weights is the steady-state
out-of-belt mean total activity.

Batches and resume
------------------
Pairs are ordered by decreasing weight and split in batches (at most
--batch-size pairs and --batch-decays decays).  Each batch is one Geant4
process driven by <outdir>/batch_NNNN.mac, with output
<outdir>/batch_NNNN_t<thread>.csv (hits) and <outdir>/batch_NNNN_runs.csv
(Volume, Isotope, NDecays of each run), and log <outdir>/batch_NNNN.log.
A batch is complete when its runs file lists all its pairs with the right
number of decays and its macro is unchanged; complete batches are skipped, the
others are deleted and rerun (unless --no-rerun).

Inputs
------
  activities.pkl, active_isotopes.pkl   from compute_activities.py
  SPENVIS file                          orbit flux (time step, belt)

Output (in --outdir, default result_postact/)
------
  allocation.csv        Volume, Isotope, Weight_Bq, NDecays, Batch
  batch_NNNN.{mac,log}, batch_NNNN_t*.csv, batch_NNNN_runs.csv
  postact_info.json     options, budget, git commit, executable, data versions
  work/                 run directory: links to the GDML files of the mass model

Usage (CLI)
-----------
    python run_postactivation.py output/activities.pkl AP8MIN.AP8.output_mean_flux_550km_SSO.txt
    python run_postactivation.py ... --dry-run        # allocation and macros only
    python run_postactivation.py ... --test 9         # timing on 9 pairs, then stop
    caffeinate -i python run_postactivation.py ... > postact.log 2>&1 &

Library use
-----------
    from run_postactivation import load_events
    hits, pairs = load_events("result_postact")
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

try:
    from activation_history import _parse_duration, _prepare_flux_and_norm, belt_mask, lag_weights
    from accumulate_spectra import primaries_per_step, source_beam_area
    from average_spectrum import weighted_average
    from compute_activities import load_outputs as load_activities, times_of
except ImportError as exc:
    sys.exit(f"ERROR: pipeline modules not found next to run_postactivation.py: {exc}")

DEFAULTS = {
    "executable":   str(HERE.parent / "Debug" / "cusp-postactivation"),
    "geometry_dir": str(HERE.parent / "gdml-mass-model"),
}
PRE_COMMANDS = ["/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year"]
HIT_DTYPES = {"RunID": "int32", "EventID": "int32", "ScintID": "int16",
              "Edep_keV": "float32", "t_ns": "float32"}


# ---------------------------------------------------------------------------
# Weights: steady-state out-of-belt activity of each (volume, isotope)
# ---------------------------------------------------------------------------

def compute_pair_weights(
    activities_pkl: str | Path,
    spenvis_file:   str | Path,
    energies:       list[float] | None = None,
    in_belt_only:   bool  = False,
    all_orbit:      bool  = False,
    belt_threshold: float = 0.0,
    duration_s:     float | None = None,
    R:              float | None = None,
    thetamax:       float | None = None,
) -> pd.DataFrame:
    """
    Steady-state mean activity [Bq] of every (volume, isotope) pair of
    active_isotopes.pkl: out-of-belt (default) or orbit-averaged (all_orbit),
    computed as in average_spectrum.compute_average_spectrum, which sums the
    same per-pair values weighted by the detector spectrum.

    Returns a DataFrame with columns Volume, Isotope, Weight_Bq, sorted by
    decreasing weight; the attrs hold the total and the parameters used.
    """
    activities_pkl = Path(activities_pkl)
    acts, active = load_activities(activities_pkl.parent, fmt=activities_pkl.suffix.lstrip("."))
    tau = times_of(acts)
    time_cols = [c for c in acts.columns if c.startswith("t_")]

    sim_energies = sorted(float(e) for e in acts.index.get_level_values("energy_MeV").unique())
    if energies is not None:
        unknown = sorted(set(map(float, energies)) - set(sim_energies))
        if unknown:
            raise ValueError(f"Energies {unknown} MeV are not in {activities_pkl.name} "
                             f"(simulated: {sim_energies}).")
        sim_energies = sorted(map(float, energies))

    # Primaries per SPENVIS time step in each energy band
    beam_area = source_beam_area(acts.attrs, R, thetamax)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # band edges above the SPENVIS range: known, reported by step 6
        nprim, _, _ = primaries_per_step(spenvis_file, sim_energies, beam_area,
                                         in_belt_only=in_belt_only)

    # Activity of each pair after one step [Bq], summed over energies
    e_level = acts.index.get_level_values("energy_MeV").astype(float)
    keep = e_level.isin(sim_energies)
    w_e = np.array([nprim[e] for e in e_level[keep]])
    pair_act = (pd.DataFrame(acts.loc[keep, time_cols].to_numpy(dtype=float) * w_e[:, None],
                             index=acts.index[keep])
                  .groupby(level=["volume", "isotope"]).sum())

    # Lag-weighted steady-state average, as in average_spectrum
    flux, _, F_mean, dt_s = _prepare_flux_and_norm(spenvis_file, in_belt_only)
    T = tau[-1] if duration_s is None else min(duration_s, tau[-1])
    select = np.ones(len(flux), bool) if all_orbit else ~belt_mask(flux, belt_threshold)
    if not select.any():
        raise ValueError("No out-of-belt steps in the SPENVIS file.")
    w_far = float(np.mean(flux) / F_mean)
    w = lag_weights(flux, F_mean, select, max(1, min(len(flux), int(T // dt_s))))
    avg = pd.Series(weighted_average(tau, pair_act.to_numpy().T, dt_s, w, w_far, T),
                    index=pair_act.index)

    pairs = [(v, i) for v, isos in active.items() for i in isos]
    weights = avg.reindex(pd.MultiIndex.from_tuples(pairs, names=["volume", "isotope"])).fillna(0.0)
    table = (weights.rename("Weight_Bq").reset_index()
             .rename(columns={"volume": "Volume", "isotope": "Isotope"}))
    table["Weight_Bq"] = table["Weight_Bq"].clip(lower=0.0)
    table = table.sort_values(["Weight_Bq", "Volume", "Isotope"],
                              ascending=[False, True, True], kind="stable").reset_index(drop=True)
    table.attrs.update({"total_Bq": float(table["Weight_Bq"].sum()),
                        "total_all_pairs_Bq": float(avg.sum()),
                        "timestep_s": dt_s, "T_s": float(T), "energies": sim_energies,
                        "in_belt_only": in_belt_only, "all_orbit": all_orbit,
                        "belt_threshold": belt_threshold,
                        "out_of_belt_fraction": float(select.mean())})
    return table


# ---------------------------------------------------------------------------
# Allocation of the decays
# ---------------------------------------------------------------------------

def allocate_decays(
    weights:  np.ndarray,
    budget:   float = 1e8,
    nmin:     int   = 1000,
    nmax:     int   = 100000,
    round_to: int   = 100,
) -> np.ndarray:
    """
    N_p = clip(s * weight_p, nmin, nmax) with the scale s chosen (bisection)
    so that sum N_p = budget, then rounded to multiples of round_to (and kept
    within [nmin, nmax]).  Zero-weight pairs get nmin.  If the budget is
    below n*nmin or above n*nmax, all pairs get nmin or nmax.
    """
    w = np.asarray(weights, dtype=float)
    lo_n = int(round(nmin / round_to)) * round_to or round_to
    hi_n = max(int(round(nmax / round_to)) * round_to, lo_n)
    if budget <= len(w) * lo_n:
        return np.full(len(w), lo_n, dtype=np.int64)
    if budget >= len(w) * hi_n:
        return np.full(len(w), hi_n, dtype=np.int64)

    def total(s: float) -> float:
        return float(np.clip(s * w, lo_n, hi_n).sum())

    pos = w[w > 0]
    s_lo, s_hi = lo_n / pos.max() * 0.5, hi_n / pos.min() * 2.0     # total(s_lo) < budget < total(s_hi)
    for _ in range(200):
        s = np.sqrt(s_lo * s_hi)
        if total(s) < budget:
            s_lo = s
        else:
            s_hi = s
        if s_hi / s_lo < 1 + 1e-12:
            break
    n = np.clip(np.round(np.clip(s_hi * w, lo_n, hi_n) / round_to) * round_to, lo_n, hi_n)
    return n.astype(np.int64)


def make_batches(table: pd.DataFrame, batch_size: int = 200, batch_decays: float = 1e6) -> list[list[int]]:
    """Row positions of `table` (already ordered) split in consecutive batches of at most
    batch_size pairs and batch_decays decays (a pair above the limit is alone)."""
    batches: list[list[int]] = []
    cur: list[int] = []
    ndec = 0
    for pos, n in enumerate(table["NDecays"].to_numpy()):
        if cur and (len(cur) >= batch_size or ndec + n > batch_decays):
            batches.append(cur)
            cur, ndec = [], 0
        cur.append(pos)
        ndec += int(n)
    if cur:
        batches.append(cur)
    return batches


def summarize_allocation(table: pd.DataFrame, nmin: int, nmax: int, top: int = 20) -> None:
    w, n = table["Weight_Bq"], table["NDecays"]
    at_min = n <= nmin
    tot = w.sum()
    print(f"Total steady-state activity of the pairs: {tot:.6g} Bq in {len(table)} pairs")
    print(f"Total decays: {int(n.sum()):,}  (nmin {nmin:,}: {int(at_min.sum())} pairs;"
          f" nmax {nmax:,}: {int((n >= nmax).sum())} pairs)")
    print(f"Weight fraction in pairs at nmax: {w[n >= nmax].sum() / tot:.4f}; "
          f"at nmin: {w[at_min].sum() / tot:.4f}")
    print(f"Top {top} pairs:")
    print(f"  {'Volume':<26}{'Isotope':<16}{'Weight [Bq]':>13}{'share':>9}{'NDecays':>10}")
    for r in table.head(top).itertuples():
        print(f"  {r.Volume:<26}{r.Isotope:<16}{r.Weight_Bq:>13.5g}{r.Weight_Bq / tot:>9.4f}{r.NDecays:>10d}")


# ---------------------------------------------------------------------------
# Macros and batch execution
# ---------------------------------------------------------------------------

def batch_name(i: int) -> str:
    return f"batch_{i:04d}"


def macro_text(prefix: Path, pairs: list[tuple[str, str, int]], pre_commands: list[str]) -> str:
    lines = ["# generated by run_postactivation.py", *PRE_COMMANDS, *pre_commands,
             f"/postact/output {prefix}"]
    for volume, isotope, n in pairs:
        lines += [f"/postact/isotope {isotope}", f"/postact/volume {volume}", f"/run/beamOn {n:d}"]
    return "\n".join(lines) + "\n"


def read_runs(path: Path) -> list[tuple[str, str, int]]:
    """Rows (Volume, Isotope, NDecays) of a runs file; [] if absent."""
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return [(r["Volume"], r["Isotope"], int(r["NDecays"])) for r in csv.DictReader(f)]


def batch_complete(outdir: Path, name: str, pairs: list[tuple[str, str, int]], macro: str) -> bool:
    mac = outdir / f"{name}.mac"
    if not mac.exists() or mac.read_text() != macro:
        return False
    return sorted(read_runs(outdir / f"{name}_runs.csv")) == sorted(pairs)


def prepare_workdir(workdir: Path, geometry_dir: Path) -> None:
    """Run directory with symlinks to every file of the mass model."""
    workdir.mkdir(parents=True, exist_ok=True)
    for f in geometry_dir.iterdir():
        link = workdir / f.name
        if f.is_file() and not link.exists():
            link.symlink_to(f)


def run_batch(exe: str, threads: int, outdir: Path, name: str, workdir: Path) -> tuple[int, float]:
    t0 = time.monotonic()
    with open(outdir / f"{name}.log", "w") as log:
        proc = subprocess.run([exe, "-t", str(threads), str(outdir / f"{name}.mac")],
                              cwd=workdir, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode, time.monotonic() - t0


def run_all(table: pd.DataFrame, batches: list[list[int]], outdir: Path, cfg: dict,
            dry_run: bool, rerun: bool, info: dict, info_path: Path) -> list[str]:
    """Write the macros and run the incomplete batches; returns the failed batch names."""
    pairs_all = list(zip(table["Volume"], table["Isotope"], table["NDecays"].astype(int)))
    exe, workdir = cfg["executable"], outdir / "work"
    if not dry_run:
        prepare_workdir(workdir, Path(cfg["geometry_dir"]))
    failed: list[str] = []
    done_decays, t_start = 0, time.monotonic()

    for i, rows in enumerate(batches):
        name = batch_name(i)
        pairs = [pairs_all[r] for r in rows]
        n_batch = sum(p[2] for p in pairs)
        macro = macro_text(outdir / name, pairs, cfg["pre_commands"])
        if not dry_run and batch_complete(outdir, name, pairs, macro):
            print(f"[{name}] complete ({len(pairs)} pairs), skipped")
            done_decays += n_batch
            continue

        if not dry_run:
            stale = [p for p in outdir.glob(f"{name}_t*.csv")] + [outdir / f"{name}_runs.csv"]
            if any(p.exists() for p in stale) and not rerun:
                print(f"[{name}] incomplete, not rerun (--no-rerun)")
                failed.append(name)
                continue
            for p in stale + [outdir / f"{name}.log"]:
                p.unlink(missing_ok=True)
        (outdir / f"{name}.mac").write_text(macro)
        if dry_run:
            continue

        print(f"[{name}] {len(pairs)} pairs, {n_batch:,} decays, "
              f"first {pairs[0][0]}/{pairs[0][1]} ...", flush=True)
        rc, wall = run_batch(exe, cfg["threads"], outdir, name, workdir)
        runs = read_runs(outdir / f"{name}_runs.csv")
        ok = rc == 0 and sorted(runs) == sorted(pairs)
        info["batches"][name] = {"pairs": len(pairs), "decays": n_batch, "wall_s": round(wall, 1),
                                 "returncode": rc, "complete": ok}
        info_path.write_text(json.dumps(info, indent=2) + "\n")
        if not ok:
            print(f"  FAILED (return code {rc}, {len(runs)}/{len(pairs)} runs; see {name}.log)")
            failed.append(name)
            continue
        done_decays += n_batch
        spent = time.monotonic() - t_start
        print(f"  done in {wall:.0f} s ({wall / n_batch * 1e3:.3g} ms/decay); "
              f"{done_decays:,} decays done, elapsed {spent / 3600:.2f} h", flush=True)
    return failed


# ---------------------------------------------------------------------------
# Timing test
# ---------------------------------------------------------------------------

def select_test_pairs(table: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    n pairs spanning the weight range, spread over the volume types
    (Scatterer, Absorber, other structure): in each type, ranks evenly
    spaced from the heaviest to the lightest pair (first the top, then the
    last, the median, ...).
    """
    kind = np.where(table["Volume"].str.contains("Scatterer"), "Scatterer",
                    np.where(table["Volume"].str.contains("Absorber"), "Absorber", "other"))
    groups = [np.flatnonzero(kind == k) for k in ("Scatterer", "Absorber", "other")]
    groups = [g for g in groups if len(g)]
    quota = [n // len(groups) + (j < n % len(groups)) for j in range(len(groups))]
    chosen: list[int] = []
    for g, q in zip(groups, quota):
        q = min(q, len(g))
        ranks = np.unique(np.round(np.linspace(0, len(g) - 1, q)).astype(int)) if q > 1 else np.array([0])
        chosen += list(g[ranks])
    return table.iloc[sorted(chosen)].reset_index(drop=True)


def run_test(table: pd.DataFrame, n_test: int, ndecays: int, outdir: Path, cfg: dict,
             total_decays: int) -> int:
    """Time n_test pairs (one process each, ndecays decays) and extrapolate to total_decays."""
    sel = select_test_pairs(table, n_test)
    tdir = outdir / "test"
    tdir.mkdir(parents=True, exist_ok=True)
    workdir = outdir / "work"
    prepare_workdir(workdir, Path(cfg["geometry_dir"]))

    def one(name: str, pair: tuple[str, str], n: int) -> float:
        for p in tdir.glob(f"{name}*"):
            p.unlink()
        (tdir / f"{name}.mac").write_text(macro_text(tdir / name, [(*pair, n)], cfg["pre_commands"]))
        rc, wall = run_batch(cfg["executable"], cfg["threads"], tdir, name, workdir)
        if rc != 0:
            print(f"  {name}: Geant4 exited with code {rc} (see {tdir / (name + '.log')})")
        return wall

    first = (sel["Volume"][0], sel["Isotope"][0])
    startup = one("probe", first, 10)
    print(f"Start-up (10 decays of {first[0]}/{first[1]}): {startup:.1f} s")
    rows = []
    for k, r in enumerate(sel.itertuples()):
        wall = one(f"test_{k:02d}", (r.Volume, r.Isotope), ndecays)
        runs = read_runs(tdir / f"test_{k:02d}_runs.csv")
        t_dec = max(wall - startup, 0.0) / ndecays
        hits = sum(sum(1 for _ in open(p)) - 1 for p in tdir.glob(f"test_{k:02d}_t*.csv"))
        rows.append((r.Volume, r.Isotope, r.Weight_Bq, wall, t_dec, hits, len(runs)))
        print(f"  {r.Volume:<28} {r.Isotope:<16}{r.Weight_Bq:>11.4g} Bq  {wall:7.1f} s  "
              f"{t_dec * 1e3:8.3g} ms/decay  {hits} hits")
    res = pd.DataFrame(rows, columns=["Volume", "Isotope", "Weight", "wall", "t_dec", "hits", "runs"])
    mean = res["t_dec"].mean()
    wsum = res["Weight"].sum()
    wmean = float((res["t_dec"] * res["Weight"]).sum() / wsum) if wsum > 0 else mean
    print(f"\nTime per decay, wall (threads={cfg['threads']}): mean {mean * 1e3:.3g} ms, "
          f"activity-weighted {wmean * 1e3:.3g} ms (range {res['t_dec'].min() * 1e3:.3g}"
          f" - {res['t_dec'].max() * 1e3:.3g} ms)")
    print(f"Full allocation ({total_decays:,} decays): {mean * total_decays / 3600:.3g} h "
          f"(mean), {wmean * total_decays / 3600:.3g} h (weighted)")
    return 0


# ---------------------------------------------------------------------------
# Run info
# ---------------------------------------------------------------------------

def git_info() -> dict:
    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True,
                              timeout=20).stdout.strip()
    try:
        return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
                "dirty": bool(git("status", "--porcelain", "--untracked-files=no"))}
    except (OSError, subprocess.SubprocessError):
        return {}


def data_versions() -> dict:
    d = os.environ.get("GEANT4_DATA_DIR") or os.environ.get("G4RADIOACTIVEDATA", "")
    if not d:
        return {}
    return {"GEANT4_DATA_DIR": d,
            "datasets": sorted(p.name for p in Path(d).glob("*") if re.search(r"\d", p.name))}


# ---------------------------------------------------------------------------
# Loading the results
# ---------------------------------------------------------------------------

def load_events(outdir: str | Path, prefix: str = "batch_") -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    All hits of a run_postactivation.py output directory.

    Returns
    -------
    hits  : DataFrame, one row per scintillator hit: Batch (int16), RunID,
            EventID (int32), DecayID (int64, unique per simulated decay:
            Batch<<48 | RunID<<32 | EventID), ScintID (int16), Edep_keV,
            t_ns (float32), Volume and Isotope (categorical).  RunID restarts
            at 0 in each batch, hence the join on (Batch, RunID).
    pairs : DataFrame indexed by (Volume, Isotope) with NDecays, the total
            number of simulated decays of the pair (decays without a hit are
            not in `hits`, so the rate per decay is hits / NDecays).
    """
    outdir = Path(outdir)
    runs_files = sorted(outdir.glob(f"{prefix}*_runs.csv"))
    if not runs_files:
        raise FileNotFoundError(f"No {prefix}*_runs.csv in {outdir}")
    runs = []
    for f in runs_files:
        b = int(re.search(rf"{re.escape(prefix)}(\d+)_runs\.csv$", f.name).group(1))
        r = pd.read_csv(f, dtype={"RunID": "int32", "NDecays": "int64"})
        r.insert(0, "Batch", b)
        runs.append(r)
    runs = pd.concat(runs, ignore_index=True)
    runs["Volume"] = runs["Volume"].astype("category")
    runs["Isotope"] = runs["Isotope"].astype("category")

    parts = []
    for f in sorted(outdir.glob(f"{prefix}*_t*.csv")):
        m = re.search(rf"{re.escape(prefix)}(\d+)_t\d+\.csv$", f.name)
        if not m:
            continue
        h = pd.read_csv(f, dtype=HIT_DTYPES)
        h.insert(0, "Batch", np.int16(int(m.group(1))))
        parts.append(h)
    if parts:
        hits = pd.concat(parts, ignore_index=True)
    else:
        hits = pd.DataFrame({c: pd.Series(dtype=t) for c, t in {"Batch": "int16", **HIT_DTYPES}.items()})
    hits["Batch"] = hits["Batch"].astype("int16")
    hits["DecayID"] = ((hits["Batch"].to_numpy(np.int64) << 48) | (hits["RunID"].to_numpy(np.int64) << 32)
                       | hits["EventID"].to_numpy(np.int64))
    keys = runs[["Batch", "RunID", "Volume", "Isotope"]].astype({"Batch": "int16"})
    hits = hits.merge(keys, on=["Batch", "RunID"], how="left", validate="many_to_one")
    if hits["Volume"].isna().any():
        warnings.warn("Some hits have no matching run in the runs files (incomplete batch?).")
    pairs = runs.groupby(["Volume", "Isotope"], observed=True)["NDecays"].sum().to_frame()
    return hits, pairs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Step 4: run cusp-postactivation for every (volume, isotope) pair.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("activities", help="activities.pkl from compute_activities.py "
                                      "(active_isotopes.pkl in the same directory)")
    p.add_argument("spenvis_file", help="SPENVIS AP9/AE9 output text file")
    p.add_argument("--outdir", default="result_postact", metavar="DIR")
    p.add_argument("--executable", default=DEFAULTS["executable"], help="cusp-postactivation executable")
    p.add_argument("--geometry-dir", dest="geometry_dir", default=DEFAULTS["geometry_dir"],
                   help="directory with the GDML files, linked into the run directory")
    p.add_argument("--threads", "-t", type=int, default=0, help="worker threads, passed as -t (0 = all cores)")
    p.add_argument("--pre-command", dest="pre_commands", action="append", default=[], metavar="CMD",
                   help="extra macro command before /postact/output (repeatable)")
    g = p.add_argument_group("weights (as in activation_history.py / average_spectrum.py)")
    g.add_argument("--energies", nargs="+", type=float, default=None, metavar="E")
    g.add_argument("--in-belt-only", action="store_true",
                   help="mean flux over the in-belt steps only (as in accumulate_spectra.py)")
    g.add_argument("--all-orbit", action="store_true", help="average over all steps instead of the out-of-belt ones")
    g.add_argument("--belt-threshold", type=float, default=0.0, metavar="FLUX",
                   help="flux above which a step is in the belt [p/cm2/s]")
    g.add_argument("--duration", default=None, metavar="T", help="mission duration, e.g. 3y (default: end of the time grid)")
    g.add_argument("--R", type=float, default=None, metavar="CM")
    g.add_argument("--thetamax", type=float, default=None, metavar="DEG")
    g = p.add_argument_group("allocation and batches")
    g.add_argument("--budget", type=float, default=1e8, help="total number of decays (default 1e8)")
    g.add_argument("--nmin", type=int, default=1000, help="minimum decays per pair")
    g.add_argument("--nmax", type=int, default=100000, help="maximum decays per pair")
    g.add_argument("--round-to", type=int, default=100, help="decays are rounded to multiples of this")
    g.add_argument("--batch-size", type=int, default=200, help="maximum pairs per Geant4 process")
    g.add_argument("--batch-decays", type=float, default=1e6, help="maximum decays per Geant4 process")
    g = p.add_argument_group("running")
    g.add_argument("--dry-run", action="store_true", help="write the allocation and the macros, do not run")
    g.add_argument("--test", type=int, default=None, metavar="N",
                   help="time N pairs spanning the weight range, extrapolate to the budget, and stop")
    g.add_argument("--test-decays", type=float, default=1e4, help="decays per pair in --test")
    g.add_argument("--no-rerun", dest="rerun", action="store_false",
                   help="do not rerun incomplete batches (report them as failed)")
    g.add_argument("--overwrite", action="store_true", help="delete the previous results of --outdir and restart")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    outdir = Path(args.outdir).resolve()
    cfg = {"executable": args.executable, "geometry_dir": args.geometry_dir,
           "threads": args.threads, "pre_commands": args.pre_commands}
    exe, geo = Path(cfg["executable"]), Path(cfg["geometry_dir"])
    if not args.dry_run and not os.access(exe, os.X_OK):
        sys.exit(f"ERROR: executable not found or not executable: {exe}")
    if not args.dry_run and (not geo.is_dir() or not any(geo.glob("*.gdml"))):
        sys.exit(f"ERROR: no GDML file in geometry directory {geo}")

    if args.overwrite and outdir.exists() and not args.dry_run:
        for p in outdir.glob("batch_*"):
            p.unlink()
        shutil.rmtree(outdir / "test", ignore_errors=True)
    outdir.mkdir(parents=True, exist_ok=True)

    table = compute_pair_weights(
        args.activities, args.spenvis_file, energies=args.energies, in_belt_only=args.in_belt_only,
        all_orbit=args.all_orbit, belt_threshold=args.belt_threshold,
        duration_s=_parse_duration(args.duration) if args.duration else None,
        R=args.R, thetamax=args.thetamax)
    table["NDecays"] = allocate_decays(table["Weight_Bq"].to_numpy(), args.budget, args.nmin,
                                       args.nmax, args.round_to)
    batches = make_batches(table, args.batch_size, args.batch_decays)
    table["Batch"] = 0
    for i, rows in enumerate(batches):
        table.loc[rows, "Batch"] = i
    table.to_csv(outdir / "allocation.csv", index=False,
                 columns=["Volume", "Isotope", "Weight_Bq", "NDecays", "Batch"], float_format="%.8g")

    a = table.attrs
    print(f"Time step {a['timestep_s']:.1f} s, T = {a['T_s']:.3g} s, "
          + ("average over all steps" if args.all_orbit
             else f"out-of-belt steps ({a['out_of_belt_fraction']:.1%} of the orbit)"))
    summarize_allocation(table, args.nmin, args.nmax)
    print(f"{len(batches)} batches (<= {args.batch_size} pairs, <= {args.batch_decays:.3g} decays)")

    info_path = outdir / "postact_info.json"
    info = {"created": dt.datetime.now().isoformat(timespec="seconds"),
            "options": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
            "executable": str(exe), "geometry_dir": str(geo), "git": git_info(),
            "geant4_data": data_versions(), "budget": args.budget,
            "total_decays": int(table["NDecays"].sum()), "n_pairs": len(table), "n_batches": len(batches),
            "weights": {k: v for k, v in a.items()}, "batches": {}}
    if info_path.exists():
        try:
            info["batches"] = json.loads(info_path.read_text()).get("batches", {})
        except json.JSONDecodeError:
            pass
    info_path.write_text(json.dumps(info, indent=2) + "\n")

    if args.test is not None:
        if args.dry_run:
            sel = select_test_pairs(table, args.test)
            print("\nTest pairs (dry run):")
            print(sel.to_string(index=False))
            return 0
        return run_test(table, args.test, int(args.test_decays), outdir, cfg, int(table["NDecays"].sum()))

    failed = run_all(table, batches, outdir, cfg, args.dry_run, args.rerun, info, info_path)
    if args.dry_run:
        print(f"\nDry run: allocation.csv and {len(batches)} macros written to {outdir}")
        return 0
    if failed:
        print(f"\nFAILED batches: {failed}")
        return 1
    print(f"\nAll {len(batches)} batches complete; load them with run_postactivation.load_events('{outdir}')")
    return 0


if __name__ == "__main__":
    sys.exit(main())
