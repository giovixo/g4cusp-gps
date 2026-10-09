"""
Parser for Geant4 activation score files.

Each CSV file corresponds to one primary energy and contains one row per
activation event: the isotope created, the volume it was created in, and
the Geant4 event ID (one primary can produce multiple activations in
different volumes, so EventID is not unique).

Output: pandas DataFrame with MultiIndex (energy_MeV, volume, isotope)
        and columns
          count           number of nuclides recorded
          efficiency      count / n_prim  [nuclides/primary]
          efficiency_err  sqrt(count) / n_prim  (Poisson)
          half_life_s     T½ in seconds, from the Lifetime column (T½ = τ·ln 2)
        Isotope names are canonical (see nuclides.py): 'Al26[228.305]'.

        df.attrs carries the run metadata needed downstream:
          'nprim'    {energy_MeV: primaries}  for every parsed file, including
                     energies with no recorded nuclide (absent from the index)
          'source'   source geometry from run_info.json (0_run.py), if found

Number of primaries: read from <dir>/run_info.json written by 0_run.py; --nprim
is only needed for score files without it, and must then match the runs.

Half-life classification used throughout:
  short    T½ <  1 000 s
  medium   1 000 s  ≤  T½  <  1 000 000 s
  long     T½  ≥  1 000 000 s

Derived accessors are provided as module-level functions rather than
subclassing DataFrame, to keep the structure transparent.

CLI USAGE:
# single file (energy inferred from name, primaries from run_info.json)
python3 activation_parser.py --file output/scorefile_10_MeV_1e+08.csv

# whole directory, save figure and pickle
python3 activation_parser.py --dir output/ --save-fig summary.png --save-pkl results.pkl

# score files without run_info.json
python3 activation_parser.py --dir ./results/ --nprim 100000000

# show top-5 instead of top-3
python3 activation_parser.py --dir ./results/ --top 5

"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from myUtilities import prettifyPlot
prettifyPlot()

from nuclides import canonical_name

_LN2 = math.log(2)

# Half-life bin edges (seconds) and labels
_HL_BINS   = [0, 1_000, 1_000_000, np.inf]
_HL_LABELS = ["Short-lived ($<$1 ks)", "Medium-lived (1 ks$–$1 Ms)", "Long-lived ($>$1 Ms)"]

_INDEX_NAMES = ["energy_MeV", "volume", "isotope"]
_COLUMNS     = ["count", "efficiency", "efficiency_err", "half_life_s"]
MANIFEST     = "run_info.json"
SCORE_GLOB   = "scorefile_*_MeV*.csv"


def _empty_frame() -> pd.DataFrame:
    idx = pd.MultiIndex.from_arrays([[], [], []], names=_INDEX_NAMES)
    return pd.DataFrame(index=idx, columns=_COLUMNS, dtype=float)


# ---------------------------------------------------------------------------
# Core parser
# ---------------------------------------------------------------------------

def parse_score_file(path: str | Path, energy_MeV: float, n_prim: int) -> pd.DataFrame:
    """
    Parse one Geant4 activation score CSV and return a tidy efficiency table.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.
    energy_MeV : float
        Primary proton energy used in this simulation run.
    n_prim : int
        Number of primaries used in this simulation run.

    Returns
    -------
    pd.DataFrame
        MultiIndex (energy_MeV, volume, isotope).
        Columns: count, efficiency, efficiency_err, half_life_s.
        Only (volume, isotope) pairs with at least one activation are stored
        (sparse — absent rows mean zero counts).  A file with no data rows
        gives an empty frame.
    """
    path = Path(path)
    raw = pd.read_csv(path, usecols=["Volume", "Isotope", "Lifetime"])
    if raw.empty:
        return _empty_frame()

    # Canonical names: Geant4 prints excitation energies with 5 decimals
    # ('Al26[228.30500]'); all pipeline steps use 3 ('Al26[228.305]').
    raw["Isotope"] = raw["Isotope"].map(canonical_name)

    # Half-life per isotope (τ is constant per isotope; take first occurrence)
    hl_map = (
        raw.groupby("Isotope")["Lifetime"]
        .first()
        .mul(_LN2)
        .rename("half_life_s")
    )

    counts = raw.groupby(["Volume", "Isotope"], sort=True).size()

    result = pd.DataFrame({"count": counts.astype(float)})
    result["efficiency"]     = result["count"] / n_prim
    result["efficiency_err"] = np.sqrt(result["count"]) / n_prim
    result["half_life_s"]    = result.index.get_level_values("Isotope").map(hl_map)

    result.index = pd.MultiIndex.from_tuples(
        [(energy_MeV, vol, iso) for vol, iso in result.index],
        names=_INDEX_NAMES,
    )

    return result


def parse_score_files(file_specs: list[tuple]) -> pd.DataFrame:
    """
    Parse multiple score files and concatenate into a single DataFrame.

    Parameters
    ----------
    file_specs : list of (path, energy_MeV, n_prim) tuples

    Returns
    -------
    pd.DataFrame
        MultiIndex (energy_MeV, volume, isotope), columns efficiency + half_life_s.
        Sorted by energy, then volume, then isotope.

    Example
    -------
    >>> specs = [
    ...     ("scorefile_10_MeV_1e8.csv",  10.0, 100_000_000),
    ...     ("scorefile_20_MeV_1e8.csv",  20.0, 100_000_000),
    ... ]
    >>> df = parse_score_files(specs)
    """
    frames = [parse_score_file(p, e, n) for p, e, n in file_specs]
    df = pd.concat(frames).sort_index()
    df.attrs["nprim"] = {float(e): int(n) for _, e, n in file_specs}
    return df


# ---------------------------------------------------------------------------
# Directory scanner
# ---------------------------------------------------------------------------

def read_manifest(directory: str | Path) -> dict | None:
    """Return the run_info.json written by 0_run.py in *directory*, or None."""
    path = Path(directory) / MANIFEST
    return json.loads(path.read_text()) if path.exists() else None


def nprim_for_file(path: Path, manifest: dict | None, n_prim: int | None) -> int:
    """Primaries for one score file: from the manifest if listed, else n_prim."""
    if manifest:
        for run in manifest.get("runs", {}).values():
            if run.get("file") == path.name:
                if n_prim is not None and int(n_prim) != int(run["nprim"]):
                    raise ValueError(
                        f"--nprim {n_prim} contradicts {MANIFEST} "
                        f"({run['nprim']} primaries for {path.name})"
                    )
                return int(run["nprim"])
    if n_prim is None:
        raise ValueError(
            f"Number of primaries for {path.name} unknown: no {MANIFEST} entry. "
            f"Pass --nprim."
        )
    return int(n_prim)


def scan_directory(
    directory: str | Path,
    n_prim: int | None = None,
    pattern: str = r"(\d+(?:\.\d+)?)_MeV",
    glob: str = SCORE_GLOB,
) -> pd.DataFrame:
    """
    Parse all score CSVs in *directory* and return the combined DataFrame.

    Parameters
    ----------
    directory : folder with the score CSVs (and run_info.json from 0_run.py)
    n_prim    : primaries per file, used only for files not listed in
                run_info.json (must agree with it otherwise)
    pattern   : regex with one capturing group that yields the energy in MeV
    glob      : filename filter

    Returns
    -------
    Combined DataFrame as from parse_score_files(), sorted by energy, with
    df.attrs['source'] set from run_info.json when available.

    Raises
    ------
    FileNotFoundError  if no matching files are found.
    ValueError         if a file has no energy tag, or n_prim is unknown.
    """
    directory = Path(directory)
    files = sorted(directory.glob(glob))
    if not files:
        raise FileNotFoundError(f"No files matching '{glob}' in {directory}")
    manifest = read_manifest(directory)

    specs = []
    for f in files:
        m = re.search(pattern, f.name)
        if m is None:
            raise ValueError(
                f"Cannot extract energy from filename '{f.name}'. "
                f"Pattern used: {pattern!r}"
            )
        specs.append((f, float(m.group(1)), nprim_for_file(f, manifest, n_prim)))
    specs.sort(key=lambda s: s[1])

    print(f"Found {len(specs)} file(s) in {directory}"
          f"{' (primaries from ' + MANIFEST + ')' if manifest else ''}:")
    for p, e, n in specs:
        print(f"  {p.name}  →  {e:g} MeV, {n:.3g} primaries")

    df = parse_score_files(specs)
    if manifest and "source" in manifest:
        df.attrs["source"] = manifest["source"]
    return df


# ---------------------------------------------------------------------------
# Half-life classification helper
# ---------------------------------------------------------------------------

def _classify_hl(half_life_s: pd.Series) -> pd.Categorical:
    """Bin a Series of half-life values into the named categories."""
    return pd.cut(half_life_s, bins=_HL_BINS, labels=_HL_LABELS, right=False)


# ---------------------------------------------------------------------------
# Helper: infer energy from filename
# ---------------------------------------------------------------------------

def energy_from_filename(path: str | Path) -> float | None:
    """
    Attempt to extract the proton energy in MeV from a filename of the form
    'scorefile_<N>_MeV_...csv'.  Returns None if the pattern is not found.

    >>> energy_from_filename("scorefile_10_MeV_1e_08.csv")
    10.0
    """
    match = re.search(r"(\d+(?:\.\d+)?)_MeV", Path(path).name)
    return float(match.group(1)) if match else None


# ---------------------------------------------------------------------------
# Derived queries
# ---------------------------------------------------------------------------

def active_volumes(df: pd.DataFrame, energy_MeV: float | None = None) -> list[str]:
    """Volumes with at least one activation (optionally filtered by energy)."""
    sub = df.loc[energy_MeV] if energy_MeV is not None else df
    return sorted(sub.index.get_level_values("volume").unique().tolist())


def isotopes_in_volume(
    df: pd.DataFrame,
    volume: str,
    energy_MeV: float | None = None,
) -> list[str]:
    """Isotopes produced in *volume* (optionally filtered by energy)."""
    if energy_MeV is not None:
        sub = df.loc[(energy_MeV, volume)]
    else:
        sub = df.xs(volume, level="volume")
    return sorted(sub.index.get_level_values("isotope").unique().tolist())


def all_isotopes(df: pd.DataFrame, energy_MeV: float | None = None) -> list[str]:
    """All isotopes produced anywhere in the mass model."""
    sub = df.loc[energy_MeV] if energy_MeV is not None else df
    return sorted(sub.index.get_level_values("isotope").unique().tolist())


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_summary(df: pd.DataFrame, save_path: str | Path | None = None) -> None:
    """
    Summary figure: number of distinct isotope species vs. primary proton energy.

    Primary y-axis   — scatter plot of the count of unique isotope species,
                       with one series per half-life category (short / medium /
                       long) plus a "total" series.  An isotope that activates
                       many times in many volumes still counts as 1 here.
    Secondary y-axis — number of distinct activated volumes at each energy.

    Parameters
    ----------
    df        : DataFrame as returned by parse_score_file / parse_score_files
    save_path : if given, save the figure to this path; otherwise display it
    """
    energies = sorted(df.index.get_level_values("energy_MeV").unique())

    # ---- build per-energy summary table ------------------------------------
    # Count distinct isotope *species* per half-life category, not efficiency sums.
    # An isotope that activates 1e6 times still counts as 1 here.
    rows = []
    for e in energies:
        sub = df.loc[e].copy().reset_index()
        # one row per unique isotope species (drop duplicate isotopes across volumes)
        iso_unique = sub.drop_duplicates("isotope")
        iso_unique = iso_unique.copy()
        iso_unique["hl_bin"] = _classify_hl(iso_unique["half_life_s"])
        bin_counts = iso_unique.groupby("hl_bin", observed=True)["isotope"].count()
        n_vol = sub["volume"].nunique()
        row = {"energy_MeV": e, "n_volumes": n_vol}
        for label in _HL_LABELS:
            row[label] = int(bin_counts.get(label, 0))
        row["Total"] = sum(row[label] for label in _HL_LABELS)
        rows.append(row)

    summary = pd.DataFrame(rows).set_index("energy_MeV")
    x = np.array(energies, dtype=float)

    # ---- figure ------------------------------------------------------------
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # colours and markers for each HL category + total
    series = [
        ("Short-lived ($<$1 ks)",      "blue", "o"),
        ("Medium-lived (1 ks$–$1 Ms)", "red", "s"),
        ("Long-lived ($>$1 Ms)",       "green", "^"),
        ("Total",              "black", "D"),
    ]

    for label, color, marker in series:
        lw     = 2.0 if label == "Total" else 1.2
        ms     = 7   if label == "Total" else 4
        ls     = "-" if label == "Total" else "--"
        zorder = 5   if label == "Total" else 3
        if label not in summary.columns:
            continue
        ax1.loglog(x, summary[label].values,
                 linestyle=ls, linewidth=lw,
                 marker=marker, markersize=ms,
                 color=color, label=label, zorder=zorder)

    # secondary axis: number of activated volumes
    ax2.loglog(x, summary["n_volumes"].values,
             "s-", linewidth=1.5,
             color='magenta')

    # x-axis: use actual energy values (not integer indices)
    ax1.set_xlabel("Energy [MeV]")
    ax1.set_ylabel("Number of isotopes")
    ax1.set_xlim(1,1000)
    ax1.set_ylim(1,1000)
   # ax1.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    
    ax2.set_ylabel("Number of activated volumes", color='magenta')
    #ax2.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax2.set_ylim(1,1000)


    # if len(energies) <= 12:
    #     ax1.set_xticks(x)
    #     ax1.set_xticklabels([f"{e:g}" for e in energies])

    # combine legends
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8, framealpha=0.85)

    ax1.grid()
    
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150)
        print(f"Figure saved to {save_path}")
    else:
        plt.show()

    plt.close(fig)


# ---------------------------------------------------------------------------
# Statistics report
# ---------------------------------------------------------------------------

def print_statistics(df: pd.DataFrame, top_n: int = 3) -> None:
    """
    Print a structured statistics report.

    For each energy:
      • Top-N most activated volumes (by total efficiency summed over isotopes)
      • Top-N most common isotopes (by total efficiency summed over volumes)

    Followed by a global summary across all energies.

    Parameters
    ----------
    df    : DataFrame as returned by parse_score_file / parse_score_files
    top_n : how many entries to show in each ranking (default 3)
    """
    energies = sorted(df.index.get_level_values("energy_MeV").unique())
    empty = sorted(set(df.attrs.get("nprim", {})) - set(energies))
    if empty:
        print(f"\nNo nuclide recorded at: {', '.join(f'{e:g}' for e in empty)} MeV")
    sep  = "─" * 60
    sep2 = "═" * 60

    print(f"\n{sep2}")
    print("  ACTIVATION STATISTICS")
    print(sep2)

    # ---- per-energy --------------------------------------------------------
    for e in energies:
        sub = df.loc[e]
        total_eff = sub["efficiency"].sum()
        n_vol     = sub.index.get_level_values("volume").nunique()
        n_iso     = sub.index.get_level_values("isotope").nunique()

        # half-life classification
        sub_reset = sub.reset_index()
        sub_reset["hl_bin"] = _classify_hl(sub_reset["half_life_s"])
        bin_totals = sub_reset.groupby("hl_bin", observed=True)["efficiency"].sum()

        print(f"\n{'─'*60}")
        print(f"  Energy: {e:g} MeV")
        print(f"{'─'*60}")
        total_err = np.sqrt(sub["count"].sum()) / sub["count"].sum() * total_eff
        print(f"  Total efficiency      : {total_eff:.4e} ± {total_err:.1e}  (activations/primary)")
        print(f"  Active volumes        : {n_vol}")
        print(f"  Distinct isotopes     : {n_iso}")
        print(f"  By half-life category :")
        for label in _HL_LABELS:
            val = bin_totals.get(label, 0.0)
            print(f"    {label:<22s}: {val:.4e}")

        # top-N volumes
        vol_eff = (
            sub.groupby(level="volume")["efficiency"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
        )
        print(f"\n  Top-{top_n} most activated volumes:")
        for rank, (vol, eff) in enumerate(vol_eff.items(), 1):
            n_iso_vol = sub.xs(vol, level="volume").index.get_level_values("isotope").nunique()
            print(f"    {rank}. {vol:<28s}  eff={eff:.4e}  ({n_iso_vol} isotope(s))")

        # top-N isotopes
        iso_eff = (
            sub.groupby(level="isotope")["efficiency"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
        )
        print(f"\n  Top-{top_n} most common isotopes:")
        hl_map = sub.reset_index().drop_duplicates("isotope").set_index("isotope")["half_life_s"]
        for rank, (iso, eff) in enumerate(iso_eff.items(), 1):
            hl  = hl_map.get(iso, float("nan"))
            cat = _classify_hl(pd.Series([hl])).iloc[0]
            print(f"    {rank}. {iso:<24s}  eff={eff:.4e}  T½={hl:.3e} s  [{cat}]")

    # ---- global summary ----------------------------------------------------
    print(f"\n{sep2}")
    print("  GLOBAL SUMMARY  (all energies combined)")
    print(sep2)

    iso_global = (
        df.groupby(level="isotope")["efficiency"]
        .sum()
        .sort_values(ascending=False)
    )
    vol_global = (
        df.groupby(level="volume")["efficiency"]
        .sum()
        .sort_values(ascending=False)
    )

    print(f"\n  Top-{top_n} volumes across all energies:")
    for rank, (vol, eff) in enumerate(vol_global.head(top_n).items(), 1):
        print(f"    {rank}. {vol:<28s}  cumulative eff={eff:.4e}")

    print(f"\n  Top-{top_n} isotopes across all energies:")
    hl_map_global = (
        df.reset_index()
        .drop_duplicates("isotope")
        .set_index("isotope")["half_life_s"]
    )
    for rank, (iso, eff) in enumerate(iso_global.head(top_n).items(), 1):
        hl  = hl_map_global.get(iso, float("nan"))
        cat = _classify_hl(pd.Series([hl])).iloc[0]
        print(f"    {rank}. {iso:<24s}  cumulative eff={eff:.4e}  T½={hl:.3e} s  [{cat}]")

    print(f"\n{sep2}\n")


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save(df: pd.DataFrame, path: str | Path) -> None:
    """
    Save the efficiency DataFrame to disk.

    Supported formats (inferred from extension):
      .pkl / .pickle  — Python pickle, fastest round-trip, Python-only
      .parquet        — columnar binary, smaller files, readable by R / Spark
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".pkl", ".pickle"):
        df.to_pickle(path)
    elif suffix == ".parquet":
        df.to_parquet(path)
    else:
        raise ValueError(
            f"Unsupported extension '{suffix}'. Use .pkl, .pickle, or .parquet."
        )


def load(path: str | Path) -> pd.DataFrame:
    """Load a DataFrame previously saved with save()."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".pkl", ".pickle"):
        return pd.read_pickle(path)
    elif suffix == ".parquet":
        return pd.read_parquet(path)
    else:
        raise ValueError(
            f"Unsupported extension '{suffix}'. Use .pkl, .pickle, or .parquet."
        )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_argparser():
    import argparse
    p = argparse.ArgumentParser(
        description="Parse Geant4 activation score files, plot, and print statistics."
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--file", metavar="CSV",
        help="Single score CSV file.",
    )
    group.add_argument(
        "--dir", metavar="DIR",
        help="Directory containing multiple score CSV files.",
    )
    p.add_argument(
        "--energy", type=float, default=None,
        help="Energy in MeV (required with --file if not inferable from filename).",
    )
    p.add_argument(
        "--nprim", type=float, default=None,
        help="Number of primaries per file; only needed for files not listed in "
             "run_info.json (default: read from run_info.json).",
    )
    p.add_argument(
        "--save-pkl", metavar="PATH",
        help="Save the combined DataFrame to a pickle file.",
    )
    p.add_argument(
        "--save-fig", metavar="PATH",
        help="Save the summary figure instead of displaying it.",
    )
    p.add_argument(
        "--top", type=int, default=3,
        help="Number of top entries in the statistics report (default: 3).",
    )
    return p


if __name__ == "__main__":
    parser = _build_argparser()
    args   = parser.parse_args()

    # ---- load data ---------------------------------------------------------
    try:
        if args.dir:
            df = scan_directory(args.dir, n_prim=args.nprim)
        else:
            energy = args.energy or energy_from_filename(args.file)
            if energy is None:
                parser.error(
                    "Cannot infer energy from filename. Provide --energy explicitly."
                )
            path = Path(args.file)
            manifest = read_manifest(path.parent)
            n_prim = nprim_for_file(path, manifest, args.nprim)
            df = parse_score_files([(path, energy, n_prim)])
            if manifest and "source" in manifest:
                df.attrs["source"] = manifest["source"]
    except (ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")

    # ---- optional save -----------------------------------------------------
    if args.save_pkl:
        save(df, args.save_pkl)
        print(f"DataFrame saved to {args.save_pkl}")

    if df.empty:
        sys.exit("No nuclide recorded in any file: nothing to plot.")

    # ---- plot --------------------------------------------------------------
    plot_summary(df, save_path=args.save_fig)

    # ---- statistics --------------------------------------------------------
    print_statistics(df, top_n=args.top)
