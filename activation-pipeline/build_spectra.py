"""
build_spectra.py
================
Step 4b of the activation pipeline: turns the hit lists written by
run_postactivation.py (result_postact/batch_NNNN_t*.csv) into detector spectra
per decay of each (volume, isotope) pair, one file per event class ("mode"),
in the format of pair_spectra.py (spectra_<mode>.npz).

THE EVENT SELECTION IS PROVISIONAL: the CUSP thresholds, coincidence window and
energy resolution are not decided yet; the defaults below are placeholders.

Selection of the events of a decay
----------------------------------
1. Optional Gaussian smearing of every scintillator deposit (off by default):
   sigma = FWHM(E)/2.355, FWHM(E) = f * sqrt(E * e_ref), f the relative FWHM at
   e_ref (--fwhm-scat / --fwhm-abs).  Fixed seed, reproducible per batch.
2. Threshold: a scintillator triggers if Edep >= --thr-scat (scatterers,
   ScintID 0-63) or --thr-abs (absorbers, ScintID 64-95).  Deposits below
   threshold are ignored entirely.
3. Coincidence window: t_first is the earliest time among the triggered
   scintillators of the decay; triggered scintillators with t - t_first >
   --window are dropped.  Simplification: in reality they would be a separate
   event; here they are ignored.  t is the time of the first deposit in each
   scintillator (the hit list has no later ones).
4. Modes (n_s, n_a = triggered scatterers, absorbers after the steps above):
     scat_single   n_s == 1, n_a == 0      E = the scatterer deposit
     abs_single    n_s == 0, n_a == 1      E = the absorber deposit
     compton       n_s == 1, n_a == 1      E = sum of the two (no geometric cuts)
     any           n_s + n_a >= 1          E = sum of all triggered deposits
5. Binning: channels of --binwidth from 0 to --emax; events with E >= emax are
   dropped (their number is in meta["n_overflow"]).

Output (in --outdir, default the input directory)
------
  spectra_<mode>.npz     see pair_spectra.py; meta holds all the parameters

Usage (CLI)
-----------
    python build_spectra.py result_postact
    python build_spectra.py result_postact --modes compton any --thr-abs 30 --window 200
    python build_spectra.py result_postact --fwhm-abs 0.07 --jobs 4

Library use
-----------
    from build_spectra import select_events, build_spectra
"""

from __future__ import annotations

import argparse
import datetime as dt
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from pair_spectra import save_pair_spectra, spectra_path

HERE = Path(__file__).resolve().parent

MODES = ("scat_single", "abs_single", "compton", "any")
N_SCATTERERS = 64                 # ScintID < 64: plastic scatterers, else GAGG absorbers
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
HIT_DTYPES = {"RunID": "int32", "EventID": "int32", "ScintID": "int16",
              "Edep_keV": "float64", "t_ns": "float64"}

DEFAULT_PARAMS = {
    "thr_scat": 5.0, "thr_abs": 20.0, "window": 500.0,
    "fwhm_scat": 0.0, "fwhm_abs": 0.0, "e_ref": 662.0,
    "emax": 2000.0, "binwidth": 1.0, "seed": 12345,
}


# ---------------------------------------------------------------------------
# Selection (one batch)
# ---------------------------------------------------------------------------

def smear(hits: pd.DataFrame, params: dict, rng: np.random.Generator) -> np.ndarray:
    """Edep of the hits after Gaussian smearing (a copy; the input is not modified)."""
    e = hits["Edep_keV"].to_numpy(np.float64).copy()
    is_scat = hits["ScintID"].to_numpy() < N_SCATTERERS
    for mask, f in ((is_scat, params["fwhm_scat"]), (~is_scat, params["fwhm_abs"])):
        if f > 0 and mask.any():
            sigma = f * np.sqrt(e[mask] * params["e_ref"]) * FWHM_TO_SIGMA
            e[mask] = e[mask] + rng.normal(0.0, 1.0, mask.sum()) * sigma
    return e


def select_events(hits: pd.DataFrame, params: dict, modes=MODES,
                  rng: np.random.Generator | None = None) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """
    Event selection on the hits of one batch (RunID, EventID, ScintID, Edep_keV, t_ns).

    Returns {mode: (run_ids, energies_keV)}: one entry per selected event
    (decay), in the order of (RunID, EventID).  Energies are not cut at emax.
    """
    empty = (np.empty(0, np.int32), np.empty(0, np.float64))
    if len(hits) == 0:
        return {m: empty for m in modes}
    if rng is None:
        rng = np.random.default_rng(params["seed"])
    run = hits["RunID"].to_numpy(np.int64)
    ev = hits["EventID"].to_numpy(np.int64)
    sid = hits["ScintID"].to_numpy()
    t = hits["t_ns"].to_numpy(np.float64)
    e = smear(hits, params, rng)
    is_scat = sid < N_SCATTERERS
    trig = e >= np.where(is_scat, params["thr_scat"], params["thr_abs"])
    run, ev, is_scat, t, e = run[trig], ev[trig], is_scat[trig], t[trig], e[trig]
    if len(e) == 0:
        return {m: empty for m in modes}

    order = np.argsort((run << 32) | ev, kind="stable")
    run, ev, is_scat, t, e = run[order], ev[order], is_scat[order], t[order], e[order]
    key = (run << 32) | ev
    start = np.concatenate(([0], np.flatnonzero(np.diff(key)) + 1))
    # coincidence window relative to the earliest triggered scintillator of the decay
    t_first = np.minimum.reduceat(t, start)
    n_per = np.diff(np.append(start, len(key)))
    keep = (t - np.repeat(t_first, n_per)) <= params["window"]
    run, key, is_scat, e = run[keep], key[keep], is_scat[keep], e[keep]
    start = np.concatenate(([0], np.flatnonzero(np.diff(key)) + 1))
    n_s = np.add.reduceat(is_scat.astype(np.int32), start)
    n_a = np.add.reduceat((~is_scat).astype(np.int32), start)
    e_sum = np.add.reduceat(e, start)
    ev_run = run[start].astype(np.int32)

    masks = {"scat_single": (n_s == 1) & (n_a == 0),
             "abs_single": (n_s == 0) & (n_a == 1),
             "compton": (n_s == 1) & (n_a == 1),
             "any": (n_s + n_a) >= 1}
    return {m: (ev_run[masks[m]], e_sum[masks[m]]) for m in modes}


def read_batch_hits(outdir: Path, batch: str) -> pd.DataFrame:
    """All thread files of one batch (batch = 'batch_0000') as one DataFrame."""
    files = sorted(f for f in outdir.glob(f"{batch}_t*.csv") if re.fullmatch(rf"{batch}_t\d+\.csv", f.name))
    parts = [pd.read_csv(f, dtype=HIT_DTYPES) for f in files]
    if not parts:
        return pd.DataFrame({c: pd.Series(dtype=t) for c, t in HIT_DTYPES.items()})
    return pd.concat(parts, ignore_index=True)


def process_batch(args: tuple) -> tuple:
    """Worker: returns (batch, runs table, {mode: (run, channel, count, n_events, n_overflow)})."""
    outdir, batch, params, modes = args
    outdir = Path(outdir)
    runs = pd.read_csv(outdir / f"{batch}_runs.csv", dtype={"RunID": "int32", "NDecays": "int64"})
    hits = read_batch_hits(outdir, batch)
    rng = np.random.default_rng(np.random.SeedSequence([params["seed"], int(batch.split("_")[1])]))
    sel = select_events(hits, params, modes, rng)
    nch = int(round(params["emax"] / params["binwidth"]))
    out = {}
    for m, (run, e) in sel.items():
        inside = e < params["emax"]
        run, e = run[inside], e[inside]
        chan = np.minimum((e / params["binwidth"]).astype(np.int64), nch - 1)
        key, cnt = np.unique(run.astype(np.int64) * nch + chan, return_counts=True)
        out[m] = (key // nch, key % nch, cnt, len(inside), int((~inside).sum()))
    return batch, runs[["RunID", "Volume", "Isotope", "NDecays"]], out


# ---------------------------------------------------------------------------
# All batches
# ---------------------------------------------------------------------------

def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(HERE), "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def build_spectra(indir: str | Path, outdir: str | Path | None = None, modes=MODES,
                  params: dict | None = None, jobs: int = 1, verbose: bool = True) -> dict[str, Path]:
    """Reads all batches of `indir`, writes spectra_<mode>.npz in `outdir`; returns the paths."""
    indir = Path(indir)
    outdir = Path(outdir) if outdir else indir
    p = {**DEFAULT_PARAMS, **(params or {})}
    modes = tuple(modes)
    nch = int(round(p["emax"] / p["binwidth"]))
    if abs(nch * p["binwidth"] - p["emax"]) > 1e-9 * p["emax"]:
        raise ValueError("emax must be a multiple of binwidth")
    edges = np.arange(nch + 1) * p["binwidth"]
    batches = sorted(f.name[:-len("_runs.csv")] for f in indir.glob("batch_*_runs.csv"))
    if not batches:
        raise FileNotFoundError(f"No batch_*_runs.csv in {indir}")
    tasks = [(str(indir), b, p, modes) for b in batches]

    pair_index: dict[tuple[str, str], int] = {}
    ndec: list[int] = []
    rows = {m: [] for m in modes}; cols = {m: [] for m in modes}; vals = {m: [] for m in modes}
    n_events = {m: 0 for m in modes}; n_over = {m: 0 for m in modes}
    t0 = time.time()
    if jobs > 1:
        pool = mp.Pool(jobs)
        it = pool.imap(process_batch, tasks)
    else:
        pool, it = None, map(process_batch, tasks)
    for k, (batch, runs, out) in enumerate(it):
        gi = np.empty(int(runs["RunID"].max()) + 1, dtype=np.int64)
        for rid, v, iso, nd in zip(runs["RunID"], runs["Volume"], runs["Isotope"], runs["NDecays"]):
            idx = pair_index.setdefault((v, iso), len(pair_index))
            if idx == len(ndec):
                ndec.append(0)
            ndec[idx] += int(nd)
            gi[rid] = idx
        for m, (run, chan, cnt, n_sel, n_ov) in out.items():
            rows[m].append(gi[run]); cols[m].append(chan); vals[m].append(cnt)
            n_events[m] += n_sel; n_over[m] += n_ov
        if verbose:
            print(f"  {batch} done ({k + 1}/{len(tasks)}, {time.time() - t0:.0f} s)", flush=True)
    if pool is not None:
        pool.close(); pool.join()

    vols = [v for v, _ in pair_index]; isos = [i for _, i in pair_index]
    meta_common = {**p, "source_dir": str(indir.resolve()), "git_commit": git_commit(),
                   "date": dt.datetime.now().isoformat(timespec="seconds"),
                   "selection": "PROVISIONAL (CUSP thresholds, window and resolution not decided)",
                   "n_batches": len(batches)}
    paths = {}
    outdir.mkdir(parents=True, exist_ok=True)
    for m in modes:
        r = np.concatenate(rows[m]) if rows[m] else np.empty(0, np.int64)
        c = np.concatenate(cols[m]) if cols[m] else np.empty(0, np.int64)
        v = np.concatenate(vals[m]) if vals[m] else np.empty(0, np.int64)
        counts = sp.coo_matrix((v, (r, c)), shape=(len(vols), nch)).tocsr()   # sums duplicates
        meta = {**meta_common, "mode": m, "n_events_selected": n_events[m] - n_over[m],
                "n_overflow": n_over[m]}
        path = spectra_path(outdir, m)
        save_pair_spectra(path, vols, isos, np.asarray(ndec, dtype=np.int64), edges, counts, meta)
        paths[m] = path
        if verbose:
            print(f"{m}: {int(counts.sum())} events in range, {n_over[m]} above emax -> {path}")
    return paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    d = DEFAULT_PARAMS
    p = argparse.ArgumentParser(description="Detector spectra per decay of each (volume, isotope) pair "
                                "from the step-4 hit lists. The selection parameters are PROVISIONAL "
                                "defaults: the CUSP values are not decided yet.")
    p.add_argument("indir", help="run_postactivation.py output directory (batch_*_runs.csv, batch_*_t*.csv)")
    p.add_argument("--outdir", default=None, help="output directory (default: indir)")
    p.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES), help="event classes to write")
    p.add_argument("--thr-scat", type=float, default=d["thr_scat"], help="scatterer threshold [keV] (provisional)")
    p.add_argument("--thr-abs", type=float, default=d["thr_abs"], help="absorber threshold [keV] (provisional)")
    p.add_argument("--window", type=float, default=d["window"],
                   help="coincidence window [ns] from the earliest triggered scintillator (provisional)")
    p.add_argument("--fwhm-scat", type=float, default=d["fwhm_scat"],
                   help="relative FWHM of the scatterers at --e-ref, scaling as sqrt(E); 0 = no smearing")
    p.add_argument("--fwhm-abs", type=float, default=d["fwhm_abs"], help="same for the absorbers")
    p.add_argument("--e-ref", type=float, default=d["e_ref"], help="reference energy of the FWHM [keV]")
    p.add_argument("--emax", type=float, default=d["emax"], help="upper edge of the spectra [keV]; higher events are dropped")
    p.add_argument("--binwidth", type=float, default=d["binwidth"], help="channel width [keV]")
    p.add_argument("--seed", type=int, default=d["seed"], help="seed of the smearing")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="processes (batches in parallel)")
    return p


def main(argv: list[str] | None = None) -> int:
    a = build_argparser().parse_args(argv)
    params = {"thr_scat": a.thr_scat, "thr_abs": a.thr_abs, "window": a.window, "fwhm_scat": a.fwhm_scat,
              "fwhm_abs": a.fwhm_abs, "e_ref": a.e_ref, "emax": a.emax, "binwidth": a.binwidth, "seed": a.seed}
    t0 = time.time()
    build_spectra(a.indir, a.outdir, a.modes, params, a.jobs)
    print(f"Done in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
