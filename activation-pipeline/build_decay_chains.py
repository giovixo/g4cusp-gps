"""
build_decay_chains.py
=====================
Read a results.pkl file produced by activation_parser.py, extract every
unique isotope across all energies and volumes, and build its full decay
chain using decay_chain_builder.py.  Output files are written to a
"DecayChains/" directory (one .dat file per isotope).

Decay modes come from the Geant4 RadioactiveDecay library and half-lives from
G4ENSDFSTATE, both found through $GEANT4_DATA_DIR (source geant4.sh) unless
given explicitly.  An existing chain file is rebuilt when it was made with
different parameters or data libraries (its '# params:' line).

Usage
-----
    python build_decay_chains.py results.pkl
    python build_decay_chains.py results.pkl --db /path/to/RadioactiveDecay6.1.2
    python build_decay_chains.py results.pkl --outdir MyChains/ --prune 0.005
    python build_decay_chains.py results.pkl --overwrite

Dependencies
------------
    activation_parser.py      (load)
    decay_chain_builder.py    (DecayDatabase, build_chain, write_chain)
    Both files must be in the same directory or on PYTHONPATH.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Imports from sibling modules
# ---------------------------------------------------------------------------

try:
    from activation_parser import load
except ImportError:
    sys.exit(
        "ERROR: activation_parser.py not found. "
        "Place it in the same directory as this script."
    )

try:
    from decay_chain_builder import (
        DecayDatabase,
        build_chain,
        write_chain,
        chain_filename,
        params_line,
        read_params_line,
        parse_isotope_name,
    )
except ImportError:
    sys.exit(
        "ERROR: decay_chain_builder.py not found. "
        "Place it in the same directory as this script."
    )


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def run(
    pkl_path:        str | Path,
    db_path:         str | Path | None = None,
    ensdf_path:      str | Path | None = None,
    outdir:          str | Path = "DecayChains",
    prune_threshold: float = 0.001,
    min_halflife:    float = 1e-6,
    overwrite:       bool  = False,
) -> None:
    pkl_path = Path(pkl_path)
    outdir   = Path(outdir)

    # --- Load activation DataFrame ---
    print(f"Loading activation data from: {pkl_path}")
    df = load(pkl_path)

    energies   = sorted(df.index.get_level_values("energy_MeV").unique())
    all_isos_raw = sorted(df.index.get_level_values("isotope").unique().tolist())

    print(f"  Energies in file : {energies}")
    print(f"  Unique isotopes  : {len(all_isos_raw)}")
    for iso in all_isos_raw:
        print(f"    {iso}")
    print()

    # --- Load decay data ---
    db = DecayDatabase(db_path, ensdf_path)
    print(f"Decay library : {db._db_path}  ({db.size} nuclide files)")
    print(f"Half-lives    : {db._ensdf_path}\n")
    params = params_line(db, prune_threshold, min_halflife)

    # --- Build chains ---
    outdir.mkdir(parents=True, exist_ok=True)

    skipped  = []
    stable   = []
    written  = []
    errors   = []

    for iso in all_isos_raw:
        outpath = outdir / chain_filename(iso)

        if outpath.exists() and not overwrite and read_params_line(outpath) == params:
            print(f"  {iso:<30s} up to date — skipping.")
            skipped.append(iso)
            continue

        print(f"  {iso:<30s} ", end="", flush=True)
        try:
            chains = build_chain(
                iso, db,
                prune_threshold=prune_threshold,
                min_halflife=min_halflife,
            )
            if not chains:
                outpath.unlink(missing_ok=True)     # stale file from an older build
                level = db.get_level(parse_isotope_name(iso))
                if level is None or level.half_life < 0 or not level.daughters:
                    msg = "stable or no decay data — no chain written."
                else:
                    msg = "all chains pruned — no chain written."
                print(msg)
                stable.append(iso)
                continue

            write_chain(iso, chains, outdir, params)
            print(f"{len(chains)} chain(s) → {outpath.name}")
            written.append(iso)

        except Exception as exc:
            print(f"ERROR: {exc}")
            errors.append((iso, str(exc)))

    # --- Summary ---
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Written  : {len(written)}")
    print(f"  Skipped  : {len(skipped)}  \t(up to date)")
    print(f"  Stable   : {len(stable)}   \t(no decay chain within the input parameters)")
    print(f"  Errors   : {len(errors)}")
    print(f"  Decay-library half-lives replaced by ENSDFSTATE (>1% off): {db.n_hl_corrected}")
    if errors:
        print()
        print("  Failed isotopes:")
        for iso, msg in errors:
            print(f"    {iso}: {msg}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build decay chains for all isotopes in a results.pkl file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "pkl",
        metavar="results.pkl",
        help="Path to the .pkl (or .parquet) file produced by activation_parser.py.",
    )
    parser.add_argument(
        "--db", default=None, metavar="DIR",
        help="Geant4 RadioactiveDecay directory (default: $G4RADIOACTIVEDATA, "
             "then $GEANT4_DATA_DIR/RadioactiveDecay*, then ./).",
    )
    parser.add_argument(
        "--ensdf", default=None, metavar="PATH",
        help="G4ENSDFSTATE directory or ENSDFSTATE.dat (default: $G4ENSDFSTATEDATA, "
             "then $GEANT4_DATA_DIR/G4ENSDFSTATE*).",
    )
    parser.add_argument(
        "--outdir", default="DecayChains", metavar="DIR",
        help="Directory where .dat chain files are written.",
    )
    parser.add_argument(
        "--prune", type=float, default=0.001, metavar="THRESHOLD",
        help="Discard chains with total branching ratio below this value.",
    )
    parser.add_argument(
        "--min-halflife", type=float, default=1e-6, metavar="SECONDS",
        help="Skip decay steps whose parent half-life is below this value.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild all files, even those up to date.",
    )
    args = parser.parse_args(argv)

    run(
        pkl_path=args.pkl,
        db_path=args.db,
        ensdf_path=args.ensdf,
        outdir=args.outdir,
        prune_threshold=args.prune,
        min_halflife=args.min_halflife,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
