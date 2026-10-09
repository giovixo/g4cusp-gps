"""
accumulate_spectra.py
=====================
For every time step in the SPENVIS orbit, compute the total background
spectrum in the detector by accumulating contributions from all isotopes,
volumes and simulation energies whose activity exceeds a given threshold.

Pipeline
--------
  activities.pkl      (from compute_activities.py)
      MultiIndex(energy_MeV, volume, isotope) × time columns
      Values: normalised activity [Bq/primary]

  spenvis flux file   (from spenvis_parser.py)
      Used to derive n_primaries(energy_band, timestep) via
      band_fluxes_to_primaries()

  result_spectra/{volume}_{isotope}_S-mode.dat
      Background spectrum [counts/s/keV/decay], one value per energy channel

Output
------
  spectra.pkl / .parquet
      dict-like: MultiIndex(time_s) × en_s channels  (float32 to save memory)
      Total background spectrum at each time step [counts/s/keV].

  count_rate.dat
      Two-column ASCII: time_s  total_count_rate [counts/s]

  count_rate.pdf / .png
      Plot of total count rate vs time.

Usage (CLI)
-----------
    python accumulate_spectra.py activities.pkl AP9MEAN.txt
    python accumulate_spectra.py activities.pkl AP9MEAN.txt \\
        --spectra-dir result_spectra/ \\
        --threshold 1.0 \\
        --R 5000 --thetamax 0.8021 --timestep 60 \\
        --outdir output/ --save-plot count_rate.pdf

Usage (library)
---------------
    from accumulate_spectra import accumulate

    spectra_df, count_rate = accumulate(
        "activities.pkl", "AP9MEAN.txt",
        spectra_dir="result_spectra/",
    )
"""

from __future__ import annotations

import argparse
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Local imports
# ---------------------------------------------------------------------------

try:
    from myUtilities import prettifyPlot   # optional personal plot style
    prettifyPlot()
except ImportError:
    pass

try:
    from compute_activities import load_outputs as load_activities
except ImportError:
    sys.exit("ERROR: compute_activities.py not found.")

try:
    from spenvis_parser import parse_spenvis, band_fluxes_to_primaries
except ImportError:
    sys.exit("ERROR: spenvis_parser.py not found.")


# ---------------------------------------------------------------------------
# Detector energy axis
# ---------------------------------------------------------------------------

def make_energy_axis() -> np.ndarray:
    """
    S-mode PI channel energy axis.

    binning_S = np.arange(18, 4114, 2)   # 2048 edges
    en_s      = binning_S[:-1] + np.diff(binning_S / 2)
    """
    binning_S = np.arange(18.0, 4112.0 + 2.0, 2.0)
    return binning_S[:-1] + np.diff(binning_S / 2.0)   # 2047 channels


EN_S = make_energy_axis()
N_CHANNELS = len(EN_S)   # 2047


# ---------------------------------------------------------------------------
# Spectrum file loader (cached)
# ---------------------------------------------------------------------------

_spectrum_cache: dict[Path, np.ndarray] = {}


def load_spectrum(path: Path) -> np.ndarray | None:
    """
    Load a background spectrum file.  Returns array of shape (N_CHANNELS,)
    in counts/s/keV/decay, or None if the file does not exist.
    """
    if path in _spectrum_cache:
        return _spectrum_cache[path]
    if not path.exists():
        return None
    spec = np.loadtxt(path)
    if len(spec) != N_CHANNELS:
        warnings.warn(
            f"Spectrum file {path.name} has {len(spec)} rows, "
            f"expected {N_CHANNELS}. Skipping."
        )
        return None
    _spectrum_cache[path] = spec
    return spec


def spectrum_path(
    spectra_dir: Path,
    volume:      str,
    isotope:     str,
    mode:        str = "S-mode",
) -> Path:
    return spectra_dir / f"{volume}_{isotope}_{mode}.dat"


# ---------------------------------------------------------------------------
# Core accumulation
# ---------------------------------------------------------------------------

def accumulate(
    activities_pkl:     str | Path,
    spenvis_file:       str | Path,
    spectra_dir:        str | Path  = "result_spectra",
    energies:           list[float] | None = None,
    activity_threshold: float       = 1.0,
    timestep:           float       = 60.0,
    R:                  float       = 5000.0,
    thetamax:           float       = 0.8021,
    in_belt_only:       bool        = False,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Accumulate background spectra over all energies, volumes, and isotopes.

    The output time grid is taken from the activities DataFrame (the decay
    time axis computed by compute_activities.py, e.g. [0.1, 1, 10, ..., 1e8] s
    after the end of irradiation).  The SPENVIS file is used solely to derive
    the orbit-mean n_primaries per energy band — a single scalar per band that
    converts normalised activities [Bq/primary] into actual activities [Bq].

    For each decay time t and each (energy_band, volume, isotope):
      1. actual_activity(t) = normalised_activity(energy_MeV, t) * n_primaries(band)
      2. if actual_activity(t) > activity_threshold:
             spectrum(t) += background_spectrum(volume, isotope) * actual_activity(t)

    The energy band index maps to simulation energies as follows: band[i]
    corresponds to simulation energy = energies[i] (lower edge convention).
    The last energy in the list (energies[-1]) has no corresponding band and
    is therefore not used.

    Parameters
    ----------
    activities_pkl     : path to activities.pkl from compute_activities.py
    spenvis_file       : path to SPENVIS AP9/AE9 output file (used for mean
                         orbit flux only; does not define the output time axis)
    spectra_dir        : directory containing {volume}_{isotope}_S-mode.dat files
    energies           : band edge energies [MeV]; must match the simulation
                         energies used in activities.pkl.  Default matches
                         the standard grid [7, 10, 15, ..., 400].
    activity_threshold : minimum actual activity [Bq] to include a contribution
    timestep           : duration of one irradiation time step [s]
    R                  : simulation geometry source sphere radius [cm]
    thetamax           : cone half-angle [degrees]
    in_belt_only       : if True, only in-belt time steps used for mean flux

    Returns
    -------
    spectra_df  : pd.DataFrame, index=time_s (decay time after irradiation [s]),
                  columns=detector energy channels [keV].
                  Total background spectrum [counts/s/keV] at each decay time.
    count_rate  : pd.Series, index=time_s.
                  Total count rate [counts/s] at each decay time (sum × Δε).
    """
    if energies is None:
        energies = [7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400]

    activities_pkl = Path(activities_pkl)
    spenvis_file   = Path(spenvis_file)
    spectra_dir    = Path(spectra_dir)

    # ---- Load normalised activities -----------------------------------
    print(f"Loading activities from: {activities_pkl}")
    activities_df, _ = load_activities(str(activities_pkl.parent),
                                       fmt=activities_pkl.suffix.lstrip("."))

    time_cols = [c for c in activities_df.columns if c.startswith("t_")]
    times_s   = np.array([float(c[2:]) for c in time_cols])
    n_times   = len(times_s)
    print(f"  {len(activities_df)} (energy, volume, isotope) entries")
    print(f"  {n_times} time points: {times_s[0]:.2e} – {times_s[-1]:.2e} s")

    # ---- Compute n_primaries per band from SPENVIS --------------------
    print(f"Loading SPENVIS flux from: {spenvis_file}")
    _, band_fluxes = parse_spenvis(
        spenvis_file, energies=energies, in_belt_only=in_belt_only
    )
    n_prim_per_band = band_fluxes_to_primaries(
        energies, band_fluxes,
        timestep=timestep, R=R, thetamax=thetamax,
    )
    # band[i] corresponds to simulation energy energies[i] (lower edge).
    # Build mapping: energy_MeV -> n_primaries
    sim_energies = energies[:-1]      # drop last edge (no band above it)
    e_to_nprim   = {float(e): n for e, n in zip(sim_energies, n_prim_per_band)}
    print(f"  {len(e_to_nprim)} energy bands → primaries:")
    for e, n in e_to_nprim.items():
        print(f"    {e:>6.1f} MeV: {n:.4e} primaries/timestep")

    # ---- Accumulate ---------------------------------------------------
    # Output array: shape (n_times, N_CHANNELS)
    total_spectra = np.zeros((n_times, N_CHANNELS), dtype=np.float64)

    missing_spectra: set[str] = set()
    missing_energies: set[float] = set()
    n_contributions = 0

    sim_energies_in_df = sorted(
        activities_df.index.get_level_values("energy_MeV").unique()
    )

    from tqdm import tqdm
    entries = list(activities_df.iterrows())
    bar = tqdm(entries, desc="Accumulating spectra", unit="entry", dynamic_ncols=True)

    for (energy_MeV, volume, isotope), row in bar:
        bar.set_postfix(E=f"{energy_MeV:g}", iso=isotope, vol=volume, refresh=False)

        # n_primaries for this simulation energy
        if energy_MeV not in e_to_nprim:
            if energy_MeV not in missing_energies:
                missing_energies.add(energy_MeV)
                bar.write(
                    f"  [WARN] No primaries mapping for energy {energy_MeV} MeV — skipped."
                )
            continue

        n_prim = e_to_nprim[energy_MeV]   # primaries per timestep (scalar)

        # Normalised activity array: shape (n_times,) [Bq/primary]
        norm_act = row[time_cols].to_numpy(dtype=float)

        # Actual activity at each time: [Bq]
        actual_act = norm_act * n_prim     # shape (n_times,)

        # Find time steps above threshold
        above = np.where(actual_act > activity_threshold)[0]
        if len(above) == 0:
            continue

        # Load background spectrum
        spec_path = spectrum_path(spectra_dir, volume, isotope)
        spec = load_spectrum(spec_path)
        if spec is None:
            key = f"{volume}/{isotope}"
            if key not in missing_spectra:
                missing_spectra.add(key)
                bar.write(f"  [WARN] No spectrum file: {spec_path.name}")
            continue

        # Accumulate: spectrum [counts/s/keV/decay] * activity [Bq=decays/s]
        # = [counts/s/keV] contribution
        # outer product: (n_above,) × (N_CHANNELS,) → (n_above, N_CHANNELS)
        total_spectra[above] += np.outer(actual_act[above], spec)
        n_contributions += len(above)

    print(f"\n  {n_contributions} (time, entry) contributions accumulated.")
    if missing_energies:
        print(f"  Energies with no band mapping: {sorted(missing_energies)}")
    if missing_spectra:
        print(f"  {len(missing_spectra)} spectrum file(s) missing.")

    # ---- Package output -----------------------------------------------
    col_names = [f"{e:.4f}" for e in EN_S]
    spectra_df = pd.DataFrame(
        total_spectra.astype(np.float32),
        index=pd.Index(times_s, name="time_s"),
        columns=col_names,
    )

    # Total count rate: sum(spectrum * dε) where dε = 2 keV (uniform bins)
    delta_e    = 2.0    # keV per channel
    count_rate = pd.Series(
        total_spectra.sum(axis=1) * delta_e,
        index=pd.Index(times_s, name="time_s"),
        name="count_rate_cps",
    )

    return spectra_df, count_rate


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_outputs(
    spectra_df:  pd.DataFrame,
    count_rate:  pd.Series,
    outdir:      str | Path = ".",
    fmt:         str        = "pkl",
) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if fmt == "parquet":
        spectra_df.to_parquet(outdir / "spectra.parquet")
    else:
        spectra_df.to_pickle(outdir / "spectra.pkl")

    np.savetxt(
        outdir / "count_rate.dat",
        np.column_stack([count_rate.index.to_numpy(), count_rate.to_numpy()]),
        header="time_s   count_rate_cps",
        fmt="%.6e",
    )
    print(f"  Spectra saved to {outdir}/spectra.{fmt}")
    print(f"  Count rate saved to {outdir}/count_rate.dat")


def load_outputs(
    outdir: str | Path = ".",
    fmt:    str        = "pkl",
) -> tuple[pd.DataFrame, pd.Series]:
    outdir = Path(outdir)
    if fmt == "parquet":
        spectra_df = pd.read_parquet(outdir / "spectra.parquet")
    else:
        spectra_df = pd.read_pickle(outdir / "spectra.pkl")

    arr = np.loadtxt(outdir / "count_rate.dat")
    count_rate = pd.Series(arr[:, 1],
                           index=pd.Index(arr[:, 0], name="time_s"),
                           name="count_rate_cps")
    return spectra_df, count_rate


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_count_rate(
    count_rate: pd.Series,
    title:      str             = "Total background count rate vs time",
    save_path:  str | Path | None = None,
) -> None:
    """
    Log-log plot of total count rate [counts/s] vs decay time [s].
    Only time points with count rate > 0 are plotted.
    """
    times = count_rate.index.to_numpy()
    rates = count_rate.to_numpy()

    mask = rates > 0
    if not mask.any():
        print("  [WARN] All count rates are zero; nothing to plot.")
        return

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.loglog(times[mask], rates[mask], "o-",
              color="steelblue", linewidth=1.3, markersize=4,
              label="total count rate")
    ax.set_xlabel("Time after irradiation (s)")
    ax.set_ylabel("Count rate (counts / s)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Count rate plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_spectra(
    spectra_df: pd.DataFrame,
    title:      str             = "Background spectra vs detector energy",
    save_path:  str | Path | None = None,
) -> None:
    """
    Log-log plot of background spectra for every time point in spectra_df.

    Each spectrum is drawn as a step curve.  Curves are coloured by decay
    time using the "magma" colormap, with the colour proportional to
    log10(time_s) so that decades are evenly spaced in colour space.
    A colorbar shows the time axis.  Time points with all-zero spectra
    are silently skipped.
    """
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    times_s  = spectra_df.index.to_numpy(dtype=float)
    energies = spectra_df.columns.to_numpy(dtype=float)

    # Only plot time points with non-zero spectra
    nonzero_mask = (spectra_df.values > 0).any(axis=1)
    if not nonzero_mask.any():
        print("  [WARN] All spectra are zero; nothing to plot.")
        return

    times_plot = times_s[nonzero_mask]
    spectra    = spectra_df.values[nonzero_mask]

    # Colour scale: log10(time) mapped linearly to [0, 1]
    log_t      = np.log10(times_plot)
    norm       = mcolors.Normalize(vmin=log_t.min(), vmax=log_t.max())
    cmap = plt.colormaps["magma"]

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, (t, spec) in enumerate(zip(times_plot, spectra)):
        color = cmap(norm(np.log10(t)))
        ax.plot(energies, np.where(spec > 0, spec, np.nan),
                drawstyle="steps-mid", linewidth=0.8,
                color=color, alpha=0.85)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Detector energy [keV]")
    ax.set_ylabel("Background rate [counts s$^{-1}$ keV$^{-1}$]")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", alpha=0.3)

    # Colorbar — use only ASCII/LaTeX-safe tick labels
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Time [s]")
    decade_ticks = np.arange(np.ceil(log_t.min()), np.floor(log_t.max()) + 1)
    cbar.set_ticks(decade_ticks)
    cbar.set_ticklabels([r"$10^{%d}$" % int(v) for v in decade_ticks])

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Spectra plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Accumulate background spectra from activities and "
                    "SPENVIS orbital flux.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("activities",   help="activities.pkl from compute_activities.py")
    p.add_argument("spenvis_file", help="SPENVIS AP9/AE9 output text file")
    p.add_argument("--spectra-dir",  default="result_spectra",  metavar="DIR",
                   help="Directory containing background spectrum .dat files.")
    p.add_argument("--energies", nargs="+", type=float,
                   default=[7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400],
                   metavar="E",
                   help="Band edge energies [MeV]. Must match simulation energies.")
    p.add_argument("--threshold", type=float, default=1.0, metavar="BQ",
                   help="Minimum actual activity [Bq] to include a contribution.")
    p.add_argument("--timestep",  type=float, default=60.0,  metavar="S")
    p.add_argument("--R",         type=float, default=5000.0, metavar="CM")
    p.add_argument("--thetamax",  type=float, default=0.8021, metavar="DEG")
    p.add_argument("--in-belt-only", action="store_true",
                   help="Use in-belt mean flux (non-zero rows only).")
    p.add_argument("--outdir",    default=".",    metavar="DIR")
    p.add_argument("--fmt",       choices=["pkl", "parquet"], default="pkl")
    p.add_argument("--save-plot", default=None,  metavar="PATH",
                   help="Save count rate plot to this path.")
    p.add_argument("--save-spectra", default=None, metavar="PATH",
                   help="Save background spectra plot to this path.")
    args = p.parse_args(argv)

    spectra_df, count_rate = accumulate(
        activities_pkl=args.activities,
        spenvis_file=args.spenvis_file,
        spectra_dir=args.spectra_dir,
        energies=args.energies,
        activity_threshold=args.threshold,
        timestep=args.timestep,
        R=args.R,
        thetamax=args.thetamax,
        in_belt_only=args.in_belt_only,
    )

    print(f"\nTotal count rate summary:")
    nonzero_cr = count_rate[count_rate > 0]
    if len(nonzero_cr):
        print(f"  Peak    : {nonzero_cr.max():.4e} counts/s  "
              f"at t={nonzero_cr.idxmax():.2e} s")
        print(f"  Mean    : {nonzero_cr.mean():.4e} counts/s  "
              f"(over {len(nonzero_cr)} active time steps)")
        print(f"  Baseline: {count_rate[count_rate == 0].shape[0]} zero time steps")
    else:
        print("  No contributions above threshold.")

    save_outputs(spectra_df, count_rate, outdir=args.outdir, fmt=args.fmt)
    plot_count_rate(count_rate, save_path=args.save_plot)
    plot_spectra(spectra_df, save_path=args.save_spectra)


if __name__ == "__main__":
    main()
