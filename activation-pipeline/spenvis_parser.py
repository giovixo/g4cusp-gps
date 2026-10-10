"""
spenvis_parser.py
=================
Parse a SPENVIS AP9/AE9 integral flux output file and return:
  a. The tabulated energy levels [MeV]
  b. The mean integral flux [protons/cm²/s] at a user-specified set of
     energies, interpolated from the tabulated values.

The mean is the time-average over all data points (including zeros when the
spacecraft is outside the belt), so it represents the orbit-averaged flux.

Interpolation is done in log–log space, consistent with the power-law shape
of trapped-proton spectra.  Energies that fall exactly on a tabulated grid
point are returned without interpolation.

Usage (CLI)
-----------
    python spenvis_parser.py AP9MEAN.txt
    python spenvis_parser.py AP9MEAN.txt --energies 7 10 15 20 30 40 50 60 70 100 150 200 300 400

Usage (library)
---------------
    from spenvis_parser import parse_spenvis

    energies, mean_fluxes = parse_spenvis("AP9MEAN.txt")
    energies, mean_fluxes = parse_spenvis("AP9MEAN.txt",
                                          energies=[7, 10, 15, 20, 30, 400])
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_spenvis(
    filepath:      str | Path,
    energies:      list[float] | None = None,
    in_belt_only:  bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a SPENVIS AP9/AE9 mean integral flux file.

    Parameters
    ----------
    filepath : path to the SPENVIS output text file
    energies : energy band edges [MeV].  The function returns the mean
               flux in each consecutive band [energies[i], energies[i+1]].
               Band flux = F(>E_lo) - F(>E_hi), where F is the integral
               flux interpolated (log–log) from the tabulated values.
               Default: [7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200,
                         300, 400]
    in_belt_only : if True, exclude time steps where all flux values are
                   zero (spacecraft outside the radiation belt) before
                   computing the mean.  The result is the average flux
                   during belt passages only.

    Returns
    -------
    energy_levels : np.ndarray
        Tabulated energy thresholds from the file header [MeV].
    band_fluxes : np.ndarray, shape (len(energies) - 1,)
        Mean flux [protons/cm²/s] in each energy band
        [energies[i], energies[i+1]].

    Raises
    ------
    ValueError : if a requested energy is above the maximum tabulated energy,
                 or if the energy level header line cannot be found.
    """
    if energies is None:
        energies = [7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400]

    filepath = Path(filepath)

    # ---- Parse header and data ------------------------------------------
    energy_levels: np.ndarray | None = None
    data_rows: list[list[float]] = []

    with open(filepath, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if "Energy levels" in line:
                    # Extract numbers after the colon
                    after_colon = line.split(":", 1)[1]
                    energy_levels = np.array(
                        [float(x) for x in re.findall(r"[\d.]+(?:e[+-]?\d+)?", after_colon)]
                    )
            else:
                stripped = line.strip()
                if not stripped:
                    continue
                values = [float(v) for v in stripped.split(",")]
                # Columns: datetime, x, y, z, flux_0, flux_1, ..., flux_n-1
                if energy_levels is not None and len(values) == 4 + len(energy_levels):
                    data_rows.append(values[4:])   # keep flux columns only

    if energy_levels is None:
        raise ValueError(
            "Could not find 'Energy levels' line in file header. "
            "Is this a SPENVIS AP9/AE9 output file?"
        )
    if not data_rows:
        raise ValueError("No data rows found in the file.")

    flux_matrix = np.array(data_rows)    # shape (n_timesteps, n_energies)

    # ---- Select rows for averaging --------------------------------------
    if in_belt_only:
        in_belt_mask = flux_matrix.any(axis=1)
        n_total    = len(flux_matrix)
        n_in_belt  = int(in_belt_mask.sum())
        if n_in_belt == 0:
            raise ValueError(
                "No in-belt time steps found (all flux rows are zero). "
                "Cannot compute in-belt mean."
            )
        flux_matrix = flux_matrix[in_belt_mask]
        print(f"In-belt filter: {n_in_belt}/{n_total} time steps retained.")

    tab_mean = flux_matrix.mean(axis=0)  # shape (n_energies,)

    # ---- Interpolate integral flux at each band edge -------------------
    # F(>E) is interpolated in log–log space (power-law spectrum).
    # Band flux = F(>E_lo) - F(>E_hi).
    energies_arr = np.asarray(energies, dtype=float)
    if len(energies_arr) < 2:
        raise ValueError("At least two energy values are required to define a band.")

    e_min = energy_levels[0]
    e_max = energy_levels[-1]

    if np.any(energies_arr > e_max):
        bad = energies_arr[energies_arr > e_max].tolist()
        raise ValueError(
            f"Requested energies {bad} MeV exceed the maximum tabulated "
            f"energy {e_max} MeV."
        )
    if np.any(energies_arr < e_min):
        import warnings
        bad = energies_arr[energies_arr < e_min].tolist()
        warnings.warn(
            f"Requested energies {bad} MeV are below the minimum tabulated "
            f"energy {e_min} MeV. F(E_min) will be used for those edges.",
            stacklevel=2,
        )
        energies_arr = np.clip(energies_arr, e_min, e_max)

    def _interp_integral_flux(e: float) -> float:
        """Log–log interpolation of F(>e) from the tabulated mean."""
        i_hi = int(np.searchsorted(energy_levels, e))
        i_lo = max(i_hi - 1, 0)
        i_hi = min(i_hi, len(energy_levels) - 1)
        if tab_mean[i_lo] == 0 and tab_mean[i_hi] == 0:
            return 0.0
        if e == energy_levels[i_lo]:          # exact tabulated point
            return float(tab_mean[i_lo])
        if tab_mean[i_lo] <= 0 or tab_mean[i_hi] <= 0:
            # One bracket is zero: linear interpolation in linear space
            t = (e - energy_levels[i_lo]) / (energy_levels[i_hi] - energy_levels[i_lo])
            return float(tab_mean[i_lo] + t * (tab_mean[i_hi] - tab_mean[i_lo]))
        log_t = (np.log(e) - np.log(energy_levels[i_lo])) / \
                (np.log(energy_levels[i_hi]) - np.log(energy_levels[i_lo]))
        return float(np.exp(np.log(tab_mean[i_lo]) +
                            log_t * (np.log(tab_mean[i_hi]) - np.log(tab_mean[i_lo]))))

    integral_at_edges = np.array([_interp_integral_flux(e) for e in energies_arr])
    band_fluxes = np.diff(integral_at_edges) * -1   # F(>E_lo) - F(>E_hi) >= 0
    band_fluxes = np.maximum(band_fluxes, 0.0)      # guard against float noise

    return energy_levels, band_fluxes



# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_differential_flux(
    energies:    list[float],
    band_fluxes: np.ndarray,
    title:       str = "Mean differential proton flux",
    save_path:   str | Path | None = None,
) -> None:
    """
    Log-log plot of mean differential flux vs energy band centre.

    Parameters
    ----------
    energies    : band edges [MeV] (length n+1)
    band_fluxes : mean flux per band [protons/cm²/s] (length n)
    title       : figure title
    save_path   : if given, save figure to this path; otherwise display
    """
    edges   = np.asarray(energies, dtype=float)
    centres = np.sqrt(edges[:-1] * edges[1:])     # geometric mean of edges
    widths  = np.diff(edges)                        # bin width [MeV]
    diff_flux = band_fluxes / widths                # [protons/cm²/s/MeV]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        centres, diff_flux,
        xerr=[centres - edges[:-1], edges[1:] - centres],
        fmt="o", capsize=4, markersize=5,
        color="steelblue", ecolor="steelblue", linewidth=1.5,
        label="AP9 mean",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Differential flux (protons / cm² / s / MeV)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Differential flux plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_total_flux_vs_time(
    filepath:  str | Path,
    title:     str = "Total integral flux vs time",
    save_path: str | Path | None = None,
) -> None:
    """
    Lin-log (linear time, log flux) plot of the total integral flux vs time.

    Time axis: seconds elapsed since the first time point (MJD → seconds).
    Flux: F(>E_min), i.e. the first (lowest-threshold) flux column.
    Zero-flux points (spacecraft outside belt) are shown at the bottom of
    the plot using a symlog y-axis so the off-belt passages are visible.

    Parameters
    ----------
    filepath  : SPENVIS output file
    title     : figure title
    save_path : if given, save figure; otherwise display
    """
    filepath = Path(filepath)

    # Re-parse raw data (MJD + flux columns) from the file
    mjd_list:       list[float] = []
    total_flux_list: list[float] = []
    n_flux_cols: int | None = None

    with open(filepath, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if "Energy levels" in line:
                    after_colon = line.split(":", 1)[1]
                    n_flux_cols = len(re.findall(r"[\d.]+(?:e[+-]?\d+)?", after_colon))
                continue
            stripped = line.strip()
            if not stripped or n_flux_cols is None:
                continue
            values = [float(v) for v in stripped.split(",")]
            if len(values) == 4 + n_flux_cols:
                mjd_list.append(values[0])
                total_flux_list.append(values[4])   # F(>E_min) = total integral flux

    if not mjd_list:
        raise ValueError("No data rows found.")

    mjd_arr   = np.array(mjd_list)
    flux_arr  = np.array(total_flux_list)
    t_s       = (mjd_arr - mjd_arr[0]) * 86400.0   # MJD → seconds, subtract t0

    fig, ax = plt.subplots(figsize=(9, 4))

    # symlog scale: linear region near zero avoids log(0) issues
    # linthresh set to 1% of the max non-zero flux
    nonzero = flux_arr[flux_arr > 0]
    linthresh = nonzero.min() * 0.5 if len(nonzero) else 1.0

    ax.semilogy(t_s, np.where(flux_arr > 0, flux_arr, np.nan),
                color="steelblue", linewidth=1.2, label="in-belt")
    # Mark the zero-flux (out-of-belt) segments at a fixed baseline
    out_mask = flux_arr == 0
    if out_mask.any():
        ax.fill_between(t_s, linthresh * 0.01, linthresh,
                        where=out_mask, color="lightgray", alpha=0.6,
                        label="out of belt")

    ax.set_xlabel("Time since start (s)")
    ax.set_ylabel("Integral flux (protons / cm² / s)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Total flux plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Primaries computation
# ---------------------------------------------------------------------------

def band_fluxes_to_primaries(
    energies:    list[float],
    band_fluxes: np.ndarray,
    timestep:    float = 60.0,
    R:           float = 5000.0,
    thetamax:    float = 0.8021,
) -> np.ndarray:
    """
    Convert mean band fluxes to number of primary protons per time step.

    The simulation geometry is a spherical source of radius R with a conical
    beam of half-angle thetamax.  The effective beam area on that sphere is
    pi * R^2 * sin^2(thetamax), so:

        N = timestep [s] * Fband [protons/cm^2/s]
            * pi * R^2 [cm^2] * sin^2(thetamax)

    Note: the formula is sometimes written without R^2 when working in
    normalised units; here R is explicit so the result is in absolute counts.

    Parameters
    ----------
    energies    : band edge energies [MeV], length n+1
    band_fluxes : mean flux per band [protons/cm^2/s], length n
                  (as returned by parse_spenvis)
    timestep    : duration of one irradiation step [s] (default: 60 s)
    R           : source sphere radius in the simulation geometry [cm]
                  (default: 5000 cm)
    thetamax    : cone half-angle [degrees] (default: 0.8021 deg)

    Returns
    -------
    n_primaries : np.ndarray, shape (n,)
        Number of primary protons per time step for each energy band.
    """
    beam_area  = np.pi * R**2 * np.sin(np.radians(thetamax))**2
    n_primaries = timestep * np.asarray(band_fluxes) * beam_area
    return n_primaries


def print_primaries(
    energies:    list[float],
    n_primaries: np.ndarray,
) -> None:
    """Print the primaries table to stdout."""
    total = n_primaries.sum()
    print(f"{'Energy band':>22s}   {'N primaries':>16s}   {'Fraction':>10s}")
    print(f"{'(MeV)':>22s}   {'(per timestep)':>16s}   {'':>10s}")
    print("-" * 56)
    for i, n in enumerate(n_primaries):
        band_str = f"{energies[i]:.4g} – {energies[i+1]:.4g}"
        frac = n / total if total > 0 else 0.0
        print(f"{band_str:>22s}   {n:>16.4e}   {frac:>10.4f}")
    print("-" * 56)
    print(f"{'Total':>22s}   {total:>16.4e}")


def plot_primaries(
    energies:    list[float],
    n_primaries: np.ndarray,
    title:       str = "Primaries per time step vs energy",
    save_path:   str | Path | None = None,
) -> None:
    """
    Bar chart of number of primaries per time step for each energy band,
    on a log y-axis with a secondary axis showing the relative fraction.

    Parameters
    ----------
    energies    : band edge energies [MeV], length n+1
    n_primaries : primaries per band (as returned by band_fluxes_to_primaries)
    title       : figure title
    save_path   : if given, save figure; otherwise display
    """
    edges    = np.asarray(energies, dtype=float)
    centres  = np.sqrt(edges[:-1] * edges[1:])
    widths   = np.diff(edges)
    total    = n_primaries.sum()
    fracs    = n_primaries / total if total > 0 else np.zeros_like(n_primaries)

    fig, ax1 = plt.subplots(figsize=(8, 5))
    bars = ax1.bar(
        centres, n_primaries,
        width=widths * 0.7,
        color="steelblue", alpha=0.8, edgecolor="white", linewidth=0.5,
        label="N primaries",
    )
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Energy (MeV)")
    ax1.set_ylabel("Primaries per time step")
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)

    ax2 = ax1.twinx()
    ax2.plot(centres, fracs, "o--", color="tomato", markersize=5,
             linewidth=1.2, label="Fraction")
    ax2.set_ylabel("Relative fraction")
    ax2.set_ylim(0, fracs.max() * 1.3)

    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=9)

    ax1.set_title(title)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Primaries plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Parse a SPENVIS AP9/AE9 flux file and report mean fluxes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", help="SPENVIS output text file.")
    parser.add_argument(
        "--energies", nargs="+", type=float, metavar="E",
        default=[7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400],
        help="Energy thresholds [MeV] at which to report the mean integral flux.",
    )
    parser.add_argument(
        "--in-belt-only", action="store_true",
        help="Average only over in-belt time steps (non-zero flux rows).",
    )
    parser.add_argument(
        "--save-diff", default=None, metavar="PATH",
        help="Save differential flux plot to this path (e.g. diff_flux.pdf)."
        " If omitted, the plot is displayed interactively.",
    )
    parser.add_argument(
        "--save-time", default=None, metavar="PATH",
        help="Save total flux vs time plot to this path."
        " If omitted, the plot is displayed interactively.",
    )
    parser.add_argument(
        "--save-prim", default=None, metavar="PATH",
        help="Save primaries-per-timestep plot to this path."
        " If omitted, the plot is displayed interactively.",
    )
    parser.add_argument(
        "--timestep", type=float, default=60.0, metavar="S",
        help="Duration of one irradiation time step [s].",
    )
    parser.add_argument(
        "--R", type=float, default=5000.0, metavar="CM",
        help="Source sphere radius in the simulation geometry [cm].",
    )
    parser.add_argument(
        "--thetamax", type=float, default=0.8021, metavar="DEG",
        help="Cone half-angle [degrees].",
    )
    args = parser.parse_args(argv)

    energy_levels, mean_fluxes = parse_spenvis(
        args.file, energies=args.energies, in_belt_only=args.in_belt_only
    )

    print(f"Tabulated energy levels (MeV): {energy_levels.tolist()}")
    print()
    avg_label = "Mean flux (in-belt)" if args.in_belt_only else "Mean band flux"
    edges = args.energies
    print(f"{'Energy band':>22s}   {avg_label:>22s}")
    print(f"{'(MeV)':>22s}   {'(protons/cm²/s)':>22s}")
    print("-" * 50)
    for i, f in enumerate(mean_fluxes):
        band_str = f"{edges[i]:.4g} – {edges[i+1]:.4g}"
        print(f"{band_str:>22s}   {f:>22.4e}")

    plot_differential_flux(
        args.energies, mean_fluxes,
        save_path=args.save_diff,
    )
    plot_total_flux_vs_time(
        args.file,
        save_path=args.save_time,
    )

    print()
    print(f"Primaries per time step  "
          f"(timestep={args.timestep:.0f} s,  R={args.R:.0f} cm,  θ_max={args.thetamax} deg)")
    print()
    n_prim = band_fluxes_to_primaries(
        args.energies, mean_fluxes,
        timestep=args.timestep, R=args.R, thetamax=args.thetamax,
    )
    print_primaries(args.energies, n_prim)
    plot_primaries(args.energies, n_prim, save_path=args.save_prim)


if __name__ == "__main__":
    main()
