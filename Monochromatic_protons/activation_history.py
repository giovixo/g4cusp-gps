"""
activation_history.py
=====================
Compute the time history of the activation-induced detector background rate.

Two modes
---------
Single-orbit mode (default)
    Convolves one pass of the SPENVIS orbital flux with the reference
    count-rate decay curve R(tau).

Long-term mode (--duration)
    Tiles the orbital flux over an arbitrary total duration under the
    assumption that the orbit repeats identically each period.  Uses a
    single FFT convolution on the full tiled array, so a 3-year run at
    60-second steps (~1.6M steps) completes in ~1-2 s.

Physical model
--------------
Each orbital timestep i is treated as an instantaneous irradiation with
fluence  phi_i = F[i] * dt  [protons/cm²].

The reference count-rate curve R(tau) was computed for a single 60-second
irradiation at the orbit-mean flux F_mean, giving n_prim = F_mean * dt * A
primaries.  By linearity of the Bateman equations:

    C[j] = sum_{k=0}^{j} R((k+1)*dt) * F[j-k] / F_mean

This is a causal convolution evaluated via FFT.

Interpolation of R(tau) onto the uniform dt grid uses lin-log (linear in
tau, logarithmic in R), which is exact for pure exponential decay.

Usage (CLI)
-----------
    # Single orbit
    python activation_history.py count_rate.dat AP9MEAN.txt

    # Long-term, 3 years, 1-week running average
    python activation_history.py count_rate.dat AP9MEAN.txt \\
        --duration 3y --avg-window 1w \\
        --save-plot history_3yr.pdf

Usage (library)
---------------
    from activation_history import compute_history, compute_long_term_history

    # Single orbit
    t, C = compute_history("count_rate.dat", "AP9MEAN.txt")

    # Long-term
    t, C, C_avg = compute_long_term_history(
        "count_rate.dat", "AP9MEAN.txt",
        duration_s=3*365.25*24*3600, avg_window_s=7*24*3600
    )
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.ndimage import uniform_filter1d
from scipy.signal import fftconvolve

try:
    from spenvis_parser import parse_spenvis
except ImportError:
    sys.exit("ERROR: spenvis_parser.py not found.")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_count_rate(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load R(tau) from a two-column count_rate.dat written by accumulate_spectra.py.

    Returns (tau [s], R [counts/s]).
    """
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def _read_total_flux(spenvis_file: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read per-timestep total integral flux (F > E_min, first flux column) and MJD
    from a SPENVIS output file.

    Returns (flux [p/cm²/s], mjd), both shape (N,).
    """
    path = Path(spenvis_file)
    n_flux_cols = None
    flux_list: list[float] = []
    mjd_list:  list[float] = []

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if "Energy levels" in line:
                    after = line.split(":", 1)[1]
                    n_flux_cols = len(re.findall(r"[\d.]+(?:e[+-]?\d+)?", after))
                continue
            stripped = line.strip()
            if not stripped or n_flux_cols is None:
                continue
            vals = [float(v) for v in stripped.split(",")]
            if len(vals) == 4 + n_flux_cols:
                mjd_list.append(vals[0])
                flux_list.append(vals[4])

    if not mjd_list:
        raise ValueError(f"No data rows found in {spenvis_file}")
    return np.array(flux_list), np.array(mjd_list)


def _parse_duration(s: str) -> float:
    """
    Parse a human-readable duration string to seconds.
    Accepted suffixes: s, m, h, d, w, mo (month=30.44d), y (year=365.25d).
    Examples: '90m', '3y', '2.5w', '1mo'
    """
    s = s.strip().lower()
    units = {
        "s":  1.0,
        "m":  60.0,
        "h":  3600.0,
        "d":  86400.0,
        "w":  7 * 86400.0,
        "mo": 30.4375 * 86400.0,
        "y":  365.25 * 86400.0,
    }
    # Try longest suffix first (so 'mo' is checked before 'm')
    for suffix in sorted(units, key=len, reverse=True):
        if s.endswith(suffix):
            return float(s[: -len(suffix)]) * units[suffix]
    raise ValueError(
        f"Cannot parse duration {s!r}. "
        f"Use a number followed by s/m/h/d/w/mo/y (e.g. '3y', '1mo', '90m')."
    )


# ---------------------------------------------------------------------------
# Interpolation
# ---------------------------------------------------------------------------

def interpolate_R_on_grid(
    tau:     np.ndarray,
    R:       np.ndarray,
    dt:      float,
    n_steps: int,
) -> np.ndarray:
    """
    Interpolate R(tau) onto a uniform grid tau_k = (k+1)*dt, k=0..n_steps-1.

    Uses lin-log interpolation (linear in tau, log in R), which is exact for
    exponential decay.  Points beyond the data range are set to 0; points
    before the first data point get R[0] (conservative nearest-neighbour).

    The result is truncated at the first step where R_grid < R_grid[0] * 1e-10
    and padded with zeros, to avoid convolving with negligibly small values.
    """
    tau_grid = (np.arange(n_steps) + 1.0) * dt

    pos = R > 0
    if not pos.any():
        return np.zeros(n_steps)

    tau_data   = tau[pos]
    log_R_data = np.log(R[pos])

    tau_query    = np.clip(tau_grid, tau_data[0], tau_data[-1])
    log_R_interp = np.interp(tau_query, tau_data, log_R_data)
    R_interp     = np.exp(log_R_interp)

    R_interp[tau_grid > tau_data[-1]] = 0.0
    R_interp[tau_grid < tau_data[0]]  = R[pos][0]

    # Truncate trailing negligible tail (speeds up FFT for long-term runs)
    threshold = R_interp[0] * 1e-12
    sig = np.where(R_interp > threshold)[0]
    if len(sig):
        R_interp[sig[-1] + 1:] = 0.0

    return R_interp


# ---------------------------------------------------------------------------
# Single-orbit computation
# ---------------------------------------------------------------------------

def _prepare_flux_and_norm(
    spenvis_file: str | Path,
    in_belt_only: bool,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Load and normalise the one-period orbital flux.

    Returns (flux_per_step, t_period, F_mean, dt).
    """
    flux_per_step, mjd = _read_total_flux(spenvis_file)
    N     = len(flux_per_step)
    dt    = float(np.median(np.diff(mjd)) * 86400.0)
    t_per = (mjd - mjd[0]) * 86400.0

    print(f"  Orbital flux: {N} timesteps, period={t_per[-1]:.1f} s, dt={dt:.1f} s")
    print(f"  In-belt steps: {(flux_per_step > 0).sum()}/{N}")

    if in_belt_only:
        ib    = flux_per_step[flux_per_step > 0]
        F_mean = ib.mean() if len(ib) else 1.0
        print(f"  F_mean (in-belt): {F_mean:.4e} p/cm²/s")
    else:
        F_mean = flux_per_step.mean()
        print(f"  F_mean (orbit-average): {F_mean:.4e} p/cm²/s")

    if F_mean == 0:
        raise ValueError("Mean flux is zero; cannot normalise.")

    return flux_per_step, t_per, F_mean, dt


def compute_history(
    count_rate_file: str | Path,
    spenvis_file:    str | Path,
    timestep:        float = 60.0,
    in_belt_only:    bool  = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the activation-induced count rate over a single orbital period.

    Parameters
    ----------
    count_rate_file : count_rate.dat from accumulate_spectra.py
    spenvis_file    : SPENVIS AP9/AE9 output file
    timestep        : orbital timestep [s] (default 60 s)
    in_belt_only    : if True, use in-belt mean flux for normalisation

    Returns
    -------
    t_orb : time axis [s], shape (N,), starting at 0
    C     : count rate [counts/s], shape (N,)
    """
    tau, R = load_count_rate(count_rate_file)
    print(f"Reference R(tau): {len(tau)} points, R_max={R.max():.4e} counts/s")

    flux_per_step, t_orb, F_mean, dt = _prepare_flux_and_norm(
        spenvis_file, in_belt_only
    )
    N = len(flux_per_step)

    R_grid    = interpolate_R_on_grid(tau, R, dt=dt, n_steps=N)
    norm_flux = flux_per_step / F_mean
    C         = np.maximum(fftconvolve(norm_flux, R_grid)[:N], 0.0)

    return t_orb, C


# ---------------------------------------------------------------------------
# Long-term computation
# ---------------------------------------------------------------------------

def compute_long_term_history(
    count_rate_file: str | Path,
    spenvis_file:    str | Path,
    duration_s:      float,
    avg_window_s:    float  = 7 * 86400.0,
    in_belt_only:    bool   = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the long-term activation count-rate history by tiling the
    one-period orbital flux.

    The orbit is assumed to repeat identically every period.  A single FFT
    convolution is performed on the full tiled time series.

    Parameters
    ----------
    count_rate_file : count_rate.dat from accumulate_spectra.py
    spenvis_file    : SPENVIS AP9/AE9 output file (defines one orbital period)
    duration_s      : total simulation duration [s]
    avg_window_s    : width of the running-average window [s] (default 1 week)
    in_belt_only    : if True, use in-belt mean flux for normalisation

    Returns
    -------
    t_total : time axis [s], shape (M,), starting at 0
    C       : instantaneous count rate [counts/s], shape (M,)
    C_avg   : running average of C over avg_window_s, shape (M,)
    """
    tau, R = load_count_rate(count_rate_file)
    print(f"Reference R(tau): {len(tau)} points, R_max={R.max():.4e} counts/s")

    flux_one_period, t_per, F_mean, dt = _prepare_flux_and_norm(
        spenvis_file, in_belt_only
    )
    N_orb    = len(flux_one_period)
    T_orb    = t_per[-1] + dt          # full period including last step
    n_periods = int(np.ceil(duration_s / T_orb))
    M         = n_periods * N_orb

    print(f"  Tiling: {n_periods} periods × {N_orb} steps = {M} total steps")
    print(f"  Total duration: {M * dt / 86400:.1f} days ({M * dt / (365.25*86400):.2f} yr)")

    # Build full time array and tiled flux
    t_total   = np.arange(M, dtype=float) * dt
    flux_full = np.tile(flux_one_period, n_periods)

    # R_grid: needs to cover at most M steps (or until R ~ 0)
    R_grid    = interpolate_R_on_grid(tau, R, dt=dt, n_steps=M)

    # FFT convolution
    norm_flux = flux_full / F_mean
    C         = np.maximum(fftconvolve(norm_flux, R_grid)[:M], 0.0)

    # Running average
    W       = max(1, int(round(avg_window_s / dt)))
    C_avg   = uniform_filter1d(C, size=W, mode="nearest")

    print(f"  Running average window: {W} steps = {W*dt/86400:.2f} days")

    return t_total, C, C_avg


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _saa_intervals(
    t: np.ndarray,
    flux: np.ndarray,
    dt: float,
) -> list[tuple[float, float]]:
    """
    Return a list of (t_start, t_end) intervals where flux > 0,
    representing SAA / radiation belt transits.
    Each interval is padded by dt/2 so it covers the full timestep bin.
    """
    in_belt = flux > 0
    intervals: list[tuple[float, float]] = []
    i = 0
    while i < len(in_belt):
        if in_belt[i]:
            j = i
            while j < len(in_belt) and in_belt[j]:
                j += 1
            intervals.append((t[i] - dt / 2.0, t[j - 1] + dt / 2.0))
            i = j
        else:
            i += 1
    return intervals


def plot_history(
    t_orb:     np.ndarray,
    C:         np.ndarray,
    flux:      np.ndarray | None = None,
    title:     str = "Activation-induced background rate vs orbital time",
    save_path: str | Path | None = None,
) -> None:
    """
    Plot the single-orbit activation count rate vs time (log y-axis).
    SAA / radiation belt transit intervals are shown as shaded regions.

    Parameters
    ----------
    flux : per-step total integral flux used to identify SAA transits.
           If None, no shading is added.
    """
    dt = float(np.median(np.diff(t_orb))) if len(t_orb) > 1 else 60.0

    fig, ax = plt.subplots(figsize=(11, 4))

    # SAA transit shading
    if flux is not None:
        for k, (t0, t1) in enumerate(_saa_intervals(t_orb, flux, dt)):
            ax.axvspan(t0, t1, color="tomato", alpha=0.18, zorder=1,
                       label="SAA transit" if k == 0 else None)

    # Count rate curve
    mask = C > 0
    if mask.any():
        ax.semilogy(t_orb[mask], C[mask], color="steelblue",
                    linewidth=1.2, label="background rate", zorder=3)

    ax.set_xlabel("Orbital time (s)")
    ax.set_ylabel("Count rate (counts / s)")
    ax.set_xlim(t_orb[0] - dt / 2.0, t_orb[-1] + dt / 2.0)
    ax.grid(True, which="both", linestyle="--", alpha=0.4, zorder=0)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title(title)
    fig.tight_layout()
    _save_or_show(fig, save_path, "History plot")


def plot_long_term_history(
    t_total:       np.ndarray,
    C:             np.ndarray,
    C_avg:         np.ndarray,
    avg_window_s:  float,
    title:         str = "Long-term activation-induced background rate",
    save_path:     str | Path | None = None,
) -> None:
    """
    Two-panel log-log plot of the long-term activation history.

    Top panel   : full instantaneous count rate (thin, semi-transparent)
                  overlaid with the running average (bold).
    Bottom panel: running average only, on a time axis in days.

    The instantaneous trace is downsampled for plotting when M > 100 000
    to keep file sizes and render times manageable.
    """
    t_days   = t_total / 86400.0
    win_days = avg_window_s / 86400.0

    # Downsample for the raw trace if very long
    MAX_PLOT = 100_000
    if len(t_total) > MAX_PLOT:
        step = len(t_total) // MAX_PLOT
        t_ds = t_days[::step]
        C_ds = C[::step]
    else:
        t_ds = t_days
        C_ds = C

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    # ---- Top panel: raw + average ----------------------------------------
    mask_ds  = C_ds > 0
    mask_avg = C_avg > 0

    if mask_ds.any():
        ax1.semilogy(t_ds[mask_ds], C_ds[mask_ds],
                     color="steelblue", linewidth=0.4, alpha=0.4,
                     label="instantaneous rate")
    if mask_avg.any():
        ax1.semilogy(t_days[mask_avg], C_avg[mask_avg],
                     color="steelblue", linewidth=1.8,
                     label=f"running avg ({win_days:.1f} d)")

    ax1.set_ylabel("Count rate (counts / s)")
    ax1.grid(True, which="both", linestyle="--", alpha=0.35)
    ax1.legend(fontsize=9, loc="upper left")
    ax1.set_title(title)

    # ---- Bottom panel: average only, linear x in days -------------------
    if mask_avg.any():
        ax2.semilogy(t_days[mask_avg], C_avg[mask_avg],
                     color="darkorange", linewidth=1.5)

    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel(f"Avg count rate (counts / s)\n[window: {win_days:.1f} d]")
    ax2.grid(True, which="both", linestyle="--", alpha=0.35)

    # Secondary x-axis in years
    ax_yr = ax2.twiny()
    ax_yr.set_xlim(np.array(ax2.get_xlim()) / 365.25)
    ax_yr.set_xlabel("Time (years)")

    fig.tight_layout()
    _save_or_show(fig, save_path, "Long-term history plot")


def _save_or_show(fig, save_path, label):
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"{label} saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Activation-induced background rate history.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("count_rate",   help="count_rate.dat from accumulate_spectra.py")
    p.add_argument("spenvis_file", help="SPENVIS AP9/AE9 output file")
    p.add_argument("--timestep",   type=float, default=60.0, metavar="S",
                   help="Orbital timestep [s].")
    p.add_argument("--in-belt-only", action="store_true",
                   help="Use in-belt mean flux for normalisation.")
    p.add_argument("--duration",   default=None, metavar="DUR",
                   help="Total simulation duration for long-term mode. "
                        "Examples: 3y, 18mo, 365d. If omitted, single orbit.")
    p.add_argument("--avg-window", default="1w", metavar="WIN",
                   help="Running-average window for long-term plot. "
                        "Examples: 1w, 30d, 1mo.")
    p.add_argument("--save-dat",   default=None, metavar="PATH",
                   help="Save time-history to a two-column ASCII file.")
    p.add_argument("--save-plot",  default=None, metavar="PATH",
                   help="Save plot to this path.")
    p.add_argument("--no-flux-overlay", action="store_true",
                   help="Omit orbital flux overlay (single-orbit mode only).")
    args = p.parse_args(argv)

    if args.duration is not None:
        # ---- Long-term mode ------------------------------------------------
        duration_s    = _parse_duration(args.duration)
        avg_window_s  = _parse_duration(args.avg_window)
        print(f"Long-term mode: duration={duration_s/86400:.1f} d, "
              f"avg_window={avg_window_s/86400:.1f} d")

        t_total, C, C_avg = compute_long_term_history(
            count_rate_file=args.count_rate,
            spenvis_file=args.spenvis_file,
            duration_s=duration_s,
            avg_window_s=avg_window_s,
            in_belt_only=args.in_belt_only,
        )

        pos = C_avg[C_avg > 0]
        if len(pos):
            print(f"\n  Peak avg rate : {pos.max():.4e} counts/s "
                  f"at t={t_total[C_avg.argmax()]/86400:.1f} d")
            print(f"  Steady-state  : {pos[-len(pos)//10:].mean():.4e} counts/s "
                  f"(last 10% of run)")

        if args.save_dat:
            np.savetxt(args.save_dat,
                       np.column_stack([t_total, C, C_avg]),
                       header="t_s   count_rate_cps   avg_count_rate_cps",
                       fmt="%.6e")
            print(f"  Time history saved to {args.save_dat}")

        plot_long_term_history(
            t_total, C, C_avg,
            avg_window_s=avg_window_s,
            save_path=args.save_plot,
        )

    else:
        # ---- Single-orbit mode ---------------------------------------------
        t_orb, C = compute_history(
            count_rate_file=args.count_rate,
            spenvis_file=args.spenvis_file,
            timestep=args.timestep,
            in_belt_only=args.in_belt_only,
        )

        pos = C[C > 0]
        if len(pos):
            print(f"\n  Peak rate : {pos.max():.4e} counts/s "
                  f"at t={t_orb[C.argmax()]:.1f} s")
            print(f"  Mean rate : {pos.mean():.4e} counts/s (active steps only)")

        if args.save_dat:
            np.savetxt(args.save_dat,
                       np.column_stack([t_orb, C]),
                       header="t_orb_s   count_rate_cps",
                       fmt="%.6e")
            print(f"  Time history saved to {args.save_dat}")

        flux_overlay = None
        if not args.no_flux_overlay:
            flux_overlay, _ = _read_total_flux(args.spenvis_file)

        plot_history(t_orb, C, flux=flux_overlay, save_path=args.save_plot)


if __name__ == "__main__":
    main()
