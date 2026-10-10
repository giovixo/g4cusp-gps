"""
compute_activities.py
=====================
Activities after an instantaneous irradiation, for every (energy, volume,
isotope), from the activation yields (results.pkl, activation_parser.py) and
the decay chains (DecayChains/, build_decay_chains.py).

Method
------
The decay chains depend only on the nuclide produced, not on where or by
which primary energy it was produced.  For every produced ("root") nuclide r
the unit activities

    a_{r,i}(t) = Σ_chains BR_c · λ_i N_i(t)       [Bq per root nucleus at t = 0]

of every nuclide i of its chains are computed once (bateman.py: analytic
solution with multiple-precision fallback, accurate to ~1e-12).  Then

    A_{E,v,i}(t) = Σ_r  eff_{E,v,r} · a_{r,i}(t)                [Bq/primary]
    σ_{E,v,i}(t) = sqrt( Σ_r (eff_err_{E,v,r} · a_{r,i}(t))² )   (Poisson, optional)

where eff is the production yield per primary of r in volume v at energy E.
The activity of every nuclide is assigned to the volume where its root was
produced.

A produced nuclide without a chain file is a nuclide with no decay data
(e.g. Ta178[0.000X]: Geant4 kills it without emission) or a missing chain
file: only its own activity, λ e^{-λt}, is computed, with the half-life of
results.pkl, and it is listed in a warning.

Reads
-----
  results.pkl            from activation_parser.py
  DecayChains/<iso>.dat  from build_decay_chains.py

Writes (in --outdir)
--------------------
  activities.pkl / .parquet
      MultiIndex(energy_MeV, volume, isotope) x time columns 't_<seconds>'.
      Normalised activity [Bq/primary].  attrs: 'times' (s), 'nprim' and
      'source' (copied from results.pkl).
  activities_err.pkl / .parquet     (with --errors)
      Same layout: 1-sigma statistical uncertainty [Bq/primary].
  unit_activities.pkl
      MultiIndex(root, isotope) x time columns: a_{r,i}(t) [Bq per root nucleus].
  active_isotopes.pkl
      {volume: set(isotope)} of the pairs whose peak activity, over all
      energies and times, is at least --threshold [Bq/primary].

Usage (CLI)
-----------
    python compute_activities.py results.pkl
    python compute_activities.py results.pkl --chains DecayChains/ --outdir output/ --errors
    python compute_activities.py results.pkl --tmin 1 --tmax 1e9 --per-decade 20
    python compute_activities.py results.pkl --times 60 3600 86400
"""

from __future__ import annotations

import argparse
import math
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from activation_parser import load as load_activation
from bateman import chain_activities
from decay_chain_builder import read_chain

_LN2 = math.log(2)


# ---------------------------------------------------------------------------
# Time grid
# ---------------------------------------------------------------------------

def default_times(tmin: float = 0.1, tmax: float = 1e9, per_decade: int = 20) -> np.ndarray:
    """Log-spaced grid from tmin to tmax [s], per_decade points per decade."""
    n = int(round(math.log10(tmax / tmin) * per_decade)) + 1
    return np.logspace(math.log10(tmin), math.log10(tmax), n)


def time_columns(times: np.ndarray) -> list[str]:
    return [f"t_{t:.6e}" for t in times]


# ---------------------------------------------------------------------------
# Unit activities per root nuclide
# ---------------------------------------------------------------------------

def unit_activities(
    root:         str,
    half_life_s:  float,
    chains_dir:   Path,
    times:        np.ndarray,
) -> tuple[dict[str, np.ndarray], bool]:
    """
    Activities [Bq per root nucleus] of all nuclides of the chains of `root`.

    Returns ({isotope: activity array}, has_chain_file).
    """
    chain_file = chains_dir / f"{root}.dat"
    if not chain_file.exists():
        lam = _LN2 / half_life_s
        return {root: lam * np.exp(-lam * times)}, False

    acc: dict[str, np.ndarray] = {}
    for br, steps in read_chain(chain_file).values():
        if not steps:
            continue
        # steps[i] = (half-life of member i, daughter of member i);
        # members: root, daughter of step 0, ..., daughter of step n-2
        # (the daughter of the last step is stable or terminal).
        members = [root] + [str(d) for _, d in steps[:-1]]
        lams = np.array([_LN2 / hl for hl, _ in steps])
        act = chain_activities(lams, times) * br
        for name, a in zip(members, act):
            if name in acc:
                acc[name] += a
            else:
                acc[name] = a.copy()
    return acc, True


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def compute_activities(
    activation_pkl:     str | Path,
    chains_dir:         str | Path = "DecayChains",
    times:              np.ndarray | None = None,
    activity_threshold: float = 1e-15,
    with_errors:        bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame, dict[str, set[str]]]:
    """
    Returns (activities, activities_err or None, unit_activities, active_isotopes).
    See the module docstring for the definitions.
    """
    times = default_times() if times is None else np.sort(np.asarray(times, dtype=float))
    chains_dir = Path(chains_dir)
    cols = time_columns(times)

    df = load_activation(activation_pkl)
    if df.empty:
        raise ValueError(f"{activation_pkl} contains no activation")

    # ---- unit activities, once per root nuclide ---------------------------
    roots = df.reset_index().drop_duplicates("isotope").set_index("isotope")["half_life_s"]
    t0 = time.monotonic()
    unit_rows: list[tuple[str, str]] = []
    unit_data: list[np.ndarray] = []
    no_chain: list[str] = []
    unit_index: dict[str, list[tuple[int, str]]] = {}
    for root, hl in roots.items():
        acts, has_chain = unit_activities(root, hl, chains_dir, times)
        if not has_chain:
            no_chain.append(root)
        unit_index[root] = []
        for iso, a in acts.items():
            unit_index[root].append((len(unit_rows), iso))
            unit_rows.append((root, iso))
            unit_data.append(a)
    U = np.array(unit_data)                                   # (n_unit, T)
    print(f"  Unit activities: {len(roots)} produced nuclides, {len(unit_rows)} "
          f"(root, nuclide) pairs, {time.monotonic() - t0:.1f} s")
    if no_chain:
        print(f"  [WARN] {len(no_chain)} nuclide(s) without chain file, own activity only: "
              f"{', '.join(sorted(no_chain))}")

    # ---- combine with the production yields (sparse product) -------------
    out_index: dict[tuple, int] = {}
    rows, cols_w, w, w2 = [], [], [], []
    has_err = "efficiency_err" in df.columns
    for (energy, volume, root), r in df.iterrows():
        eff = r["efficiency"]
        err = r["efficiency_err"] if has_err else 0.0
        for u, iso in unit_index[root]:
            o = out_index.setdefault((energy, volume, iso), len(out_index))
            rows.append(o); cols_w.append(u); w.append(eff); w2.append(err * err)
    shape = (len(out_index), len(unit_rows))
    W = sparse.csr_matrix((w, (rows, cols_w)), shape=shape)
    A = W @ U

    idx = pd.MultiIndex.from_tuples(list(out_index), names=["energy_MeV", "volume", "isotope"])
    activities = pd.DataFrame(A, index=idx, columns=cols).sort_index()
    activities_err = None
    if with_errors:
        W2 = sparse.csr_matrix((w2, (rows, cols_w)), shape=shape)
        activities_err = pd.DataFrame(np.sqrt(W2 @ (U * U)), index=idx, columns=cols).sort_index()

    meta = {"times": times.tolist(),
            "nprim": df.attrs.get("nprim", {}),
            "source": df.attrs.get("source")}
    activities.attrs.update(meta)
    if activities_err is not None:
        activities_err.attrs.update(meta)

    unit_df = pd.DataFrame(U, index=pd.MultiIndex.from_tuples(unit_rows, names=["root", "isotope"]),
                           columns=cols)
    unit_df.attrs["times"] = times.tolist()

    # ---- active isotopes ----------------------------------------------------
    peak = activities.max(axis=1)
    active_isotopes: dict[str, set[str]] = {}
    for (energy, volume, iso), p in peak.items():
        if p >= activity_threshold:
            active_isotopes.setdefault(volume, set()).add(iso)

    return activities, activities_err, unit_df, active_isotopes


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def _write(df: pd.DataFrame, path_stem: Path, fmt: str) -> Path:
    if fmt == "parquet":
        path = path_stem.with_suffix(".parquet")
        df.to_parquet(path)
    else:
        path = path_stem.with_suffix(".pkl")
        df.to_pickle(path)
    return path


def save_outputs(
    activities_df:   pd.DataFrame,
    active_isotopes: dict[str, set[str]],
    outdir:          str | Path = ".",
    fmt:             str = "pkl",
    activities_err:  pd.DataFrame | None = None,
    unit_df:         pd.DataFrame | None = None,
) -> None:
    """Save the outputs to outdir (fmt: 'pkl' or 'parquet' for the activity tables)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"  Activities saved to {_write(activities_df, outdir / 'activities', fmt)}")
    if activities_err is not None:
        print(f"  Uncertainties saved to {_write(activities_err, outdir / 'activities_err', fmt)}")
    if unit_df is not None:
        unit_df.to_pickle(outdir / "unit_activities.pkl")
        print(f"  Unit activities saved to {outdir / 'unit_activities.pkl'}")
    with open(outdir / "active_isotopes.pkl", "wb") as f:
        pickle.dump(active_isotopes, f)
    print(f"  Active isotopes saved to {outdir / 'active_isotopes.pkl'}")


def load_outputs(
    outdir: str | Path = ".",
    fmt:    str = "pkl",
) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    """Load activities and active_isotopes previously saved by save_outputs()."""
    outdir = Path(outdir)
    if fmt == "parquet":
        activities_df = pd.read_parquet(outdir / "activities.parquet")
    else:
        activities_df = pd.read_pickle(outdir / "activities.pkl")
    with open(outdir / "active_isotopes.pkl", "rb") as f:
        active_isotopes = pickle.load(f)
    return activities_df, active_isotopes


def times_of(df: pd.DataFrame) -> np.ndarray:
    """Evaluation times [s] of an activities table (from its column labels)."""
    return np.array([float(c[2:]) for c in df.columns if c.startswith("t_")])


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_summary(
    activities_df:   pd.DataFrame,
    active_isotopes: dict[str, set[str]],
    top_n:           int = 5,
) -> None:
    """Human-readable summary: entries, and the dominant nuclides at a few times."""
    times = times_of(activities_df)
    energies = sorted(activities_df.index.get_level_values("energy_MeV").unique())
    print(f"\n{'=' * 60}\nACTIVITY SUMMARY\n{'=' * 60}")
    print(f"Energies [MeV]: {[f'{e:g}' for e in energies]}")
    print(f"Entries (energy, volume, isotope): {len(activities_df)}")
    print(f"Time grid: {len(times)} points, {times[0]:.3g} – {times[-1]:.3g} s")
    print(f"Volumes with active isotopes: {len(active_isotopes)}, "
          f"(volume, isotope) pairs: {sum(len(s) for s in active_isotopes.values())}")
    by_iso = activities_df.groupby(level="isotope").sum()       # summed over E and volumes
    for t_target in (60.0, 3600.0, 86400.0, 3.15576e7):
        c = activities_df.columns[np.argmin(np.abs(np.log(times / t_target)))]
        top = by_iso[c].sort_values(ascending=False).head(top_n)
        print(f"  t = {float(c[2:]):9.3g} s: " +
              ", ".join(f"{i} {v:.2e}" for i, v in top.items()) +
              "   [Bq/primary, summed over energies]")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Activities after an instantaneous irradiation, from results.pkl "
                    "and the decay chains.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("pkl", metavar="results.pkl", help="Activation table from activation_parser.py.")
    p.add_argument("--chains", default="DecayChains", metavar="DIR",
                   help="Directory with the decay chain files.")
    p.add_argument("--times", nargs="+", type=float, default=None, metavar="T",
                   help="Explicit evaluation times [s] (overrides the grid options).")
    p.add_argument("--tmin", type=float, default=0.1, help="Grid start [s].")
    p.add_argument("--tmax", type=float, default=1e9, help="Grid end [s].")
    p.add_argument("--per-decade", type=int, default=20, help="Grid points per decade.")
    p.add_argument("--threshold", type=float, default=1e-15, metavar="BQ",
                   help="Min peak activity [Bq/primary] for active_isotopes.")
    p.add_argument("--errors", action="store_true",
                   help="Also compute and save the statistical uncertainties.")
    p.add_argument("--outdir", default=".", metavar="DIR", help="Output directory.")
    p.add_argument("--fmt", choices=["pkl", "parquet"], default="pkl",
                   help="Storage format for the activity tables.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    times = (np.asarray(args.times) if args.times
             else default_times(args.tmin, args.tmax, args.per_decade))

    print(f"Activation data : {args.pkl}")
    print(f"Decay chains    : {args.chains}")
    print(f"Time grid       : {len(times)} points, {times[0]:.3g} – {times[-1]:.3g} s")
    activities, act_err, unit_df, active = compute_activities(
        args.pkl, args.chains, times, args.threshold, args.errors)
    print_summary(activities, active)
    save_outputs(activities, active, outdir=args.outdir, fmt=args.fmt,
                 activities_err=act_err, unit_df=unit_df)


if __name__ == "__main__":
    main()
