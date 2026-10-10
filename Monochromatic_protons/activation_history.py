"""
activation_history.py
=====================
Compute the time history of the activation-induced detector background rate.

Two modes
---------
Single-period mode (default)
    Convolves one pass of the SPENVIS orbital flux file with the reference
    count-rate curve R(tau), starting from no activation.

Long-term mode (--duration)
    Tiles the SPENVIS flux over an arbitrary total duration under the
    assumption that it repeats identically each period.  Uses a single FFT
    convolution on the full tiled array, so a 3-year run at 60-second steps
    (~1.6M steps) completes in a few seconds.

Physical model
--------------
R(tau) [counts/s] (accumulate_spectra.py) is the rate at time tau after one
time step dt of irradiation at the mean flux F_mean.  The SPENVIS flux F[i] is
taken as constant during each step i, so by linearity of the Bateman
equations the rate at the end of step j is

    C[j] = sum_{k=0}^{j} K_k * F[j-k] / F_mean,
    K_k  = (1/dt) * integral_{k dt}^{(k+1) dt} R(tau) dtau,

the step-averaged kernel (it is also the rate averaged over step j if each
step is an instantaneous irradiation at its start).  Sampling R at the lags
(k+1) dt instead would miss the decays within one step of the irradiation:
under constant flux a nuclide of mean life tau_m would get only
x e^-x / (1 - e^-x), x = dt / tau_m, of its equilibrium rate.

R is log-linear between its grid points (exact for a single exponential),
constant before the first point and zero after the last.

F is the integral flux above the lowest simulation energy E_min, F(>E_min),
interpolated between the SPENVIS levels: protons below E_min are not simulated
and produce no activity, and their time profile along the orbit differs from
that of the activating protons (in the AP8MIN 550 km file the peak-to-mean
ratio of the total flux, F(>0.1 MeV), is twice that of F(>5 MeV)).
accumulate_spectra.py records E_min in count_rate.dat ('profile_emin_MeV:');
older files without it fall back to the total flux, with a warning.

Out-of-belt averages
--------------------
A step is in the belt when the total integral flux (first SPENVIS column, not
the profile F(>E_min)) is above --belt-threshold (default 0).  The running average of the long-term
mode and the summary rates use the out-of-belt steps only.

Usage (CLI)
-----------
    # Single period
    python activation_history.py count_rate.dat AP9MEAN.txt

    # Long-term, 3 years, 1-week running average
    python activation_history.py count_rate.dat AP9MEAN.txt \\
        --duration 3y --avg-window 1w \\
        --save-plot history_3yr.pdf

Usage (library)
---------------
    from activation_history import compute_history, compute_long_term_history

    # Single period
    t, C, out = compute_history("count_rate.dat", "AP9MEAN.txt")

    # Long-term
    t, C, C_avg, out = compute_long_term_history(
        "count_rate.dat", "AP9MEAN.txt",
        duration_s=3*365.25*24*3600, avg_window_s=7*24*3600
    )
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from functools import lru_cache
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


def count_rate_info(path: str | Path) -> dict:
    """
    Source and energy band of a count_rate.dat, from the 'source:', 'mode:' and
    'band_keV:' header lines written by accumulate_spectra.py when it used a
    spectra_<mode>.npz or a band, and the energy of the flux time profile
    ('profile_emin_MeV:', written since the profile is F(>E_min)).  Returns {}
    for the plain files, otherwise {'source', 'mode', 'emin', 'emax',
    'profile_emin', 'label'} (those present); label is a short description
    ('compton, 20-100 keV') for plot titles.
    """
    info: dict = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            key, _, val = line[1:].partition(":")
            key, val = key.strip(), val.strip()
            if key == "source":
                info["source"] = val
            elif key == "mode":
                info["mode"] = val
            elif key == "band_keV":
                lo, hi = val.split()
                info["emin"] = None if lo == "None" else float(lo)
                info["emax"] = None if hi == "None" else float(hi)
            elif key == "profile_emin_MeV":
                info["profile_emin"] = None if val == "None" else float(val)
    if info:
        lo, hi = info.get("emin"), info.get("emax")
        band = ("" if lo is None and hi is None else
                f"{lo:g}-{hi:g} keV" if lo is not None and hi is not None else
                f">= {lo:g} keV" if lo is not None else f"<= {hi:g} keV")
        info["label"] = ", ".join(x for x in (info.get("mode"), band) if x)
    return info


@lru_cache(maxsize=8)
def _read_flux_table(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Energy levels [MeV], integral flux F(>E) of every time step [p/cm²/s]
    (shape N x n_levels) and MJD (shape N) of a SPENVIS output file.  Cached:
    the arrays are read-only.
    """
    levels = None
    rows: list[list[float]] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if "Energy levels" in line:
                    after = line.split(":", 1)[1]
                    levels = [float(x) for x in re.findall(r"[\d.]+(?:e[+-]?\d+)?", after)]
                continue
            stripped = line.strip()
            if not stripped or levels is None:
                continue
            vals = [float(v) for v in stripped.split(",")]
            if len(vals) == 4 + len(levels):
                rows.append([vals[0]] + vals[4:])
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    table = np.array(rows)
    out = (np.array(levels), table[:, 1:], table[:, 0])
    for a in out:
        a.setflags(write=False)
    return out


def read_flux_table(spenvis_file: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(energy levels [MeV], F(>E) per time step [p/cm²/s], MJD) of a SPENVIS file."""
    return _read_flux_table(str(Path(spenvis_file).resolve()))


def _read_total_flux(spenvis_file: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read per-timestep total integral flux (F > E_min, first flux column) and MJD
    from a SPENVIS output file.  It defines the belt passages.

    Returns (flux [p/cm²/s], mjd), both shape (N,).
    """
    _, flux, mjd = read_flux_table(spenvis_file)
    return flux[:, 0], mjd


def flux_profile(spenvis_file: str | Path, emin: float | None) -> np.ndarray:
    """
    Integral flux F(>emin) [p/cm²/s] of every time step: the time profile of
    the activating protons, with emin the lowest simulation energy.  Log-log
    interpolation between the tabulated levels (linear where a bracket is
    zero), as in spenvis_parser.parse_spenvis.  emin=None, or emin at or below
    the first level, gives the first column (the total flux).
    """
    levels, flux, _ = read_flux_table(spenvis_file)
    if emin is None or emin <= levels[0]:
        return flux[:, 0]
    if emin > levels[-1]:
        raise ValueError(f"profile energy {emin:g} MeV is above the last level "
                         f"({levels[-1]:g} MeV) of {spenvis_file}")
    hi = int(np.searchsorted(levels, emin))
    if levels[hi] == emin:
        return flux[:, hi]
    lo = hi - 1
    f0, f1 = flux[:, lo], flux[:, hi]
    s = (np.log(emin) - np.log(levels[lo])) / (np.log(levels[hi]) - np.log(levels[lo]))
    lin = (emin - levels[lo]) / (levels[hi] - levels[lo])
    with np.errstate(divide="ignore", invalid="ignore"):
        loglog = np.exp(np.log(f0) + s * (np.log(f1) - np.log(f0)))
    return np.where((f0 > 0) & (f1 > 0), loglog, f0 + lin * (f1 - f0))


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
# Decay-curve integrals and kernel
# ---------------------------------------------------------------------------

def cumulative_integral(
    tau:    np.ndarray,
    values: np.ndarray,
    t:      np.ndarray | float,
) -> np.ndarray:
    """
    Integral from 0 to each t of every column of values (shape n_tau x n_col,
    or n_tau), with the curve log-linear between the tau points (exact for a
    single exponential), constant before tau[0] and zero after tau[-1].

    Returns shape (len(t), n_col).
    """
    tau = np.asarray(tau, dtype=float)
    v   = np.asarray(values, dtype=float).reshape(len(tau), -1)
    t   = np.atleast_1d(np.asarray(t, dtype=float))

    h  = np.diff(tau)[:, None]
    y0 = v[:-1]
    y1 = v[1:]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r = np.log(y1 / y0)
        loglin = (y0 > 0) & (y1 > 0) & (np.abs(r) > 1e-9)
        seg = np.where(loglin, h * (y1 - y0) / r, 0.5 * h * (y0 + y1))
    I_nodes = tau[0] * v[0] + np.vstack([np.zeros(v.shape[1]), np.cumsum(seg, axis=0)])

    out    = np.empty((len(t), v.shape[1]))
    before = t <= tau[0]
    after  = t >= tau[-1]
    mid    = ~before & ~after
    out[before] = t[before, None] * v[0]
    out[after]  = I_nodes[-1]
    if mid.any():
        n = np.searchsorted(tau, t[mid], side="right") - 1
        s = ((t[mid] - tau[n]) / (tau[n + 1] - tau[n]))[:, None]
        a, b, hh, rr = y0[n], y1[n], h[n], r[n]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            part = np.where(loglin[n], hh * a * np.expm1(s * rr) / rr,
                            hh * (a * s + 0.5 * (b - a) * s**2))
        out[mid] = I_nodes[n] + part
    return out


def step_kernel(
    tau:     np.ndarray,
    R:       np.ndarray,
    dt:      float,
    n_steps: int,
) -> np.ndarray:
    """
    Step-averaged kernel K_k = (1/dt) * integral_{k dt}^{(k+1) dt} R dtau,
    k = 0 .. n_steps-1, for every column of R.  Shape (n_steps, n_col).
    """
    I = cumulative_integral(tau, R, np.arange(n_steps + 1) * dt)
    return np.diff(I, axis=0) / dt


def belt_mask(flux: np.ndarray, threshold: float = 0.0) -> np.ndarray:
    """In-belt time steps: total integral flux above threshold [p/cm²/s]."""
    return np.asarray(flux) > threshold


def lag_weights(
    flux:   np.ndarray,
    F_mean: float,
    select: np.ndarray,
    n_lags: int,
) -> np.ndarray:
    """
    w_k = mean over the selected steps j of F[j-k] / F_mean, k < n_lags, for a
    flux that repeats with the period of the array (circular correlation).

    The average over the selected steps of the long-term rate is then
    sum_k w_k K_k.  For n_lags larger than the period the weights repeat.
    """
    f   = np.asarray(flux, dtype=float) / F_mean
    sel = np.asarray(select, dtype=float)
    N   = len(f)
    c   = np.fft.irfft(np.fft.rfft(sel) * np.conj(np.fft.rfft(f)), N) / sel.sum()
    return np.resize(c, n_lags)


# ---------------------------------------------------------------------------
# Single-orbit computation
# ---------------------------------------------------------------------------

def _prepare_flux_and_norm(
    spenvis_file: str | Path,
    in_belt_only: bool,
    profile_emin: float | None = None,
) -> tuple[np.ndarray, np.ndarray, float, float, np.ndarray]:
    """
    Load and normalise the one-period orbital flux.

    profile_emin : energy [MeV] of the time profile F(>profile_emin), the
                   lowest simulation energy; None uses the total flux (first
                   column), as before the profile energy was recorded.

    Returns (flux_per_step, t_period, F_mean, dt, total_flux): the time profile,
    its mean (over the in-belt steps if in_belt_only) and the total flux, which
    defines the belt passages.
    """
    total, mjd = _read_total_flux(spenvis_file)
    flux_per_step = flux_profile(spenvis_file, profile_emin)
    N     = len(flux_per_step)
    dt    = float(np.median(np.diff(mjd)) * 86400.0)
    t_per = (mjd - mjd[0]) * 86400.0

    print(f"  Orbital flux: {N} timesteps, period={t_per[-1]:.1f} s, dt={dt:.1f} s")
    print(f"  In-belt steps: {(total > 0).sum()}/{N}")
    print(f"  Time profile: F(>{profile_emin:g} MeV)" if profile_emin is not None
          else "  Time profile: total flux (first column)")

    if in_belt_only:
        ib    = flux_per_step[total > 0]
        F_mean = ib.mean() if len(ib) else 1.0
        print(f"  F_mean (in-belt): {F_mean:.4e} p/cm²/s")
    else:
        F_mean = flux_per_step.mean()
        print(f"  F_mean (orbit-average): {F_mean:.4e} p/cm²/s")

    if F_mean == 0:
        raise ValueError("Mean flux is zero; cannot normalise.")

    return flux_per_step, t_per, F_mean, dt, total


def profile_energy(recorded: float | None, override: float | None, source: str) -> float | None:
    """
    Energy of the flux time profile: override if given, else the one recorded
    by accumulate_spectra.py; files written before it was recorded fall back
    to the total flux, with a warning.
    """
    if override is not None:
        return override
    if recorded is None:
        warnings.warn(f"{source} does not record the energy of the flux time profile "
                      f"(written before it was): the total flux is used. Rerun "
                      f"accumulate_spectra.py, or pass the lowest simulation energy.",
                      stacklevel=3)
    return recorded


def compute_history(
    count_rate_file: str | Path,
    spenvis_file:    str | Path,
    in_belt_only:    bool  = False,
    belt_threshold:  float = 0.0,
    profile_emin:    float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the activation-induced count rate over one SPENVIS period,
    starting from no activation.

    Parameters
    ----------
    count_rate_file : count_rate.dat from accumulate_spectra.py
    spenvis_file    : SPENVIS AP9/AE9 output file
    in_belt_only    : if True, use in-belt mean flux for normalisation
    belt_threshold  : total flux above which a step is in the belt [p/cm²/s]
    profile_emin    : energy [MeV] of the flux time profile (default: the one
                      recorded in count_rate_file)

    Returns
    -------
    t   : time axis [s], shape (N,), starting at 0
    C   : count rate at the end of each step [counts/s], shape (N,)
    out : out-of-belt steps, boolean shape (N,)
    """
    tau, R = load_count_rate(count_rate_file)
    print(f"Reference R(tau): {len(tau)} points, R_max={R.max():.4e} counts/s")

    emin = profile_energy(count_rate_info(count_rate_file).get("profile_emin"),
                          profile_emin, str(count_rate_file))
    flux_per_step, t_orb, F_mean, dt, total = _prepare_flux_and_norm(
        spenvis_file, in_belt_only, emin
    )
    N = len(flux_per_step)

    K   = step_kernel(tau, R, dt, N)[:, 0]
    C   = np.maximum(fftconvolve(flux_per_step / F_mean, K)[:N], 0.0)
    out = ~belt_mask(total, belt_threshold)

    return t_orb, C, out


# ---------------------------------------------------------------------------
# Long-term computation
# ---------------------------------------------------------------------------

def compute_long_term_history(
    count_rate_file: str | Path,
    spenvis_file:    str | Path,
    duration_s:      float,
    avg_window_s:    float  = 7 * 86400.0,
    in_belt_only:    bool   = False,
    belt_threshold:  float  = 0.0,
    profile_emin:    float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    belt_threshold  : total flux above which a step is in the belt [p/cm²/s]
    profile_emin    : energy [MeV] of the flux time profile (default: the one
                      recorded in count_rate_file)

    Returns
    -------
    t_total : time axis [s], shape (M,), starting at 0
    C       : count rate at the end of each step [counts/s], shape (M,)
    C_avg   : running average of C over the out-of-belt steps within
              avg_window_s, shape (M,)
    out     : out-of-belt steps, boolean shape (M,)
    """
    tau, R = load_count_rate(count_rate_file)
    print(f"Reference R(tau): {len(tau)} points, R_max={R.max():.4e} counts/s")

    emin = profile_energy(count_rate_info(count_rate_file).get("profile_emin"),
                          profile_emin, str(count_rate_file))
    flux_one_period, t_per, F_mean, dt, total_one = _prepare_flux_and_norm(
        spenvis_file, in_belt_only, emin
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

    # FFT convolution with the step-averaged kernel
    K = step_kernel(tau, R, dt, M)[:, 0]
    C = np.maximum(fftconvolve(flux_full / F_mean, K)[:M], 0.0)

    # Running average over the out-of-belt steps
    out = ~belt_mask(np.tile(total_one, n_periods), belt_threshold)
    W   = max(1, int(round(avg_window_s / dt)))
    num = uniform_filter1d(C * out, size=W, mode="nearest")
    den = uniform_filter1d(out.astype(float), size=W, mode="nearest")
    with np.errstate(divide="ignore", invalid="ignore"):
        C_avg = np.where(den > 0, num / den, np.nan)

    print(f"  Running average window: {W} steps = {W*dt/86400:.2f} days "
          f"(out-of-belt steps: {out[:N_orb].mean():.1%})")

    return t_total, C, C_avg, out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _saa_intervals(
    t: np.ndarray,
    in_belt: np.ndarray,
    dt: float,
) -> list[tuple[float, float]]:
    """
    Return a list of (t_start, t_end) intervals of in-belt steps,
    representing SAA / radiation belt transits.
    Each interval is padded by dt/2 so it covers the full timestep bin.
    """
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
    in_belt:   np.ndarray | None = None,
    title:     str = "Activation-induced background rate vs orbital time",
    save_path: str | Path | None = None,
    label:     str = "",
) -> None:
    """
    Plot the single-orbit activation count rate vs time (log y-axis).
    SAA / radiation belt transit intervals are shown as shaded regions.

    Parameters
    ----------
    in_belt : in-belt steps (shaded).  If None, no shading is added.
    label   : spectra mode and energy band of the count rate, added to the title.
    """
    dt = float(np.median(np.diff(t_orb))) if len(t_orb) > 1 else 60.0

    fig, ax = plt.subplots(figsize=(11, 4))

    # SAA transit shading
    if in_belt is not None:
        for k, (t0, t1) in enumerate(_saa_intervals(t_orb, in_belt, dt)):
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
    ax.set_title(title + (f" [{label}]" if label else ""))
    fig.tight_layout()
    _save_or_show(fig, save_path, "History plot")


def plot_long_term_history(
    t_total:       np.ndarray,
    C:             np.ndarray,
    C_avg:         np.ndarray,
    avg_window_s:  float,
    title:         str = "Long-term activation-induced background rate",
    save_path:     str | Path | None = None,
    label:         str = "",
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
    mask_avg = np.nan_to_num(C_avg) > 0

    if mask_ds.any():
        ax1.semilogy(t_ds[mask_ds], C_ds[mask_ds],
                     color="steelblue", linewidth=0.4, alpha=0.4,
                     label="instantaneous rate")
    if mask_avg.any():
        ax1.semilogy(t_days[mask_avg], C_avg[mask_avg],
                     color="steelblue", linewidth=1.8,
                     label=f"out-of-belt running avg ({win_days:.1f} d)")

    ax1.set_ylabel("Count rate (counts / s)")
    ax1.grid(True, which="both", linestyle="--", alpha=0.35)
    ax1.legend(fontsize=9, loc="upper left")
    ax1.set_title(title + (f" [{label}]" if label else ""))

    # ---- Bottom panel: average only, linear x in days -------------------
    if mask_avg.any():
        ax2.semilogy(t_days[mask_avg], C_avg[mask_avg],
                     color="darkorange", linewidth=1.5)

    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel(f"Out-of-belt avg rate (counts / s)\n[window: {win_days:.1f} d]")
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
    p.add_argument("count_rate",   help="count_rate.dat from accumulate_spectra.py "
                                        "(its spectra source and energy band, set there "
                                        "with --spectra-file/--emin/--emax, carry over)")
    p.add_argument("spenvis_file", help="SPENVIS AP9/AE9 output file")
    p.add_argument("--in-belt-only", action="store_true",
                   help="Use in-belt mean flux for normalisation.")
    p.add_argument("--belt-threshold", type=float, default=0.0, metavar="FLUX",
                   help="Total flux [p/cm²/s] above which a step is in the belt.")
    p.add_argument("--profile-emin", type=float, default=None, metavar="MEV",
                   help="Energy of the flux time profile F(>E) [MeV]; default: the one "
                        "recorded in count_rate.dat (the lowest simulation energy).")
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

    info = count_rate_info(args.count_rate)
    if "source" in info:
        print(f"Count rate: {info['label'] or 'whole axis'} from {info.get('source')}")
    label = info.get("label", "")
    dat_note = (f"\nsource: {info['source']}\nmode: {info.get('mode')}\nband_keV: {info.get('emin')} {info.get('emax')}"
                if "source" in info else "")

    if args.duration is not None:
        # ---- Long-term mode ------------------------------------------------
        duration_s    = _parse_duration(args.duration)
        avg_window_s  = _parse_duration(args.avg_window)
        print(f"Long-term mode: duration={duration_s/86400:.1f} d, "
              f"avg_window={avg_window_s/86400:.1f} d")

        t_total, C, C_avg, out = compute_long_term_history(
            count_rate_file=args.count_rate,
            spenvis_file=args.spenvis_file,
            duration_s=duration_s,
            avg_window_s=avg_window_s,
            in_belt_only=args.in_belt_only,
            belt_threshold=args.belt_threshold,
            profile_emin=args.profile_emin,
        )

        flux_one, _ = _read_total_flux(args.spenvis_file)
        last = np.arange(len(C)) >= len(C) - len(flux_one)
        if out.any():
            print(f"\n  Peak out-of-belt avg rate : {np.nanmax(C_avg):.4e} counts/s "
                  f"at t={t_total[np.nanargmax(C_avg)]/86400:.1f} d")
            print(f"  Out-of-belt mean, last period : {C[out & last].mean():.4e} counts/s")
            print(f"  All-step mean, last period    : {C[last].mean():.4e} counts/s")

        if args.save_dat:
            np.savetxt(args.save_dat,
                       np.column_stack([t_total, C, C_avg]),
                       header="t_s   count_rate_cps   out_of_belt_avg_count_rate_cps" + dat_note,
                       fmt="%.6e")
            print(f"  Time history saved to {args.save_dat}")

        plot_long_term_history(
            t_total, C, C_avg,
            avg_window_s=avg_window_s,
            save_path=args.save_plot,
            label=label,
        )

    else:
        # ---- Single-orbit mode ---------------------------------------------
        t_orb, C, out = compute_history(
            count_rate_file=args.count_rate,
            spenvis_file=args.spenvis_file,
            in_belt_only=args.in_belt_only,
            belt_threshold=args.belt_threshold,
            profile_emin=args.profile_emin,
        )

        print(f"\n  Peak rate        : {C.max():.4e} counts/s at t={t_orb[C.argmax()]:.1f} s")
        if out.any():
            print(f"  Out-of-belt mean : {C[out].mean():.4e} counts/s")
        print(f"  All-step mean    : {C.mean():.4e} counts/s")

        if args.save_dat:
            np.savetxt(args.save_dat,
                       np.column_stack([t_orb, C]),
                       header="t_orb_s   count_rate_cps" + dat_note,
                       fmt="%.6e")
            print(f"  Time history saved to {args.save_dat}")

        plot_history(t_orb, C, in_belt=None if args.no_flux_overlay else ~out,
                     save_path=args.save_plot, label=label)


if __name__ == "__main__":
    main()
