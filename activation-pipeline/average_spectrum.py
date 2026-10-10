"""
average_spectrum.py
===================
Compute the steady-state out-of-belt background spectrum after a long
mission duration, and identify the isotopes responsible for the most
prominent spectral lines.

Physical basis
--------------
S(tau, E) [counts/s/keV] is the spectrum at decay time tau after one SPENVIS
time step dt of irradiation at the mean flux F_mean (accumulate_spectra.py).
With the SPENVIS flux F repeating periodically and constant during each step,
the spectrum at the end of step j is (as in activation_history.py; F is the
flux above the lowest simulation energy, F(>E_min), recorded in spectra.pkl,
and the belt passages are defined by the total flux)

    S_j(E) = sum_k K_k(E) F[j-k] / F_mean,
    K_k(E) = (1/dt) * integral_{k dt}^{(k+1) dt} S(tau, E) dtau.

Its average over the out-of-belt steps j, after a mission of duration T, is

    <S(E)> = sum_{k < T/dt} w_k K_k(E),   w_k = mean_{j out of belt} F[j-k] / F_mean.

w_k is computed exactly for the lags within one SPENVIS period; for longer
lags it is replaced by its period mean (the orbit-mean F / F_mean), since the
kernel then varies little within a period.  With --all-orbit the average is
over all steps and reduces to (<F>/F_mean / dt) * integral_0^T S dtau.

Nuclides with half-lives much longer than T have not reached their steady
state and count with the activity built up to T.  Without --duration, T is the
end of the decay-time grid of compute_activities.py (1e9 s by default).

Per-isotope decomposition
-------------------------
If activities_pkl and a spectra source are given (--spectra-dir with the .dat
files, or --spectra-file with a spectra_<mode>.npz from build_spectra.py; the
same source and mode as in accumulate_spectra.py), the activity of each (volume,
isotope) after one step, A_{v,i}(tau) = sum_j N_j a_{j,v,i}(tau), is averaged
in the same way and multiplied by its spectrum per decay F_{v,i}(E):

    <S_i(E)> = sum_v <A_{v,i}> F_{v,i}(E).

The isotopes sum to <S(E)> when accumulate_spectra.py ran without an activity
threshold.  The dominant isotope at each spectral peak is the one with the
largest <S_i(E_peak)>.

Usage (CLI)
-----------
    # Total spectrum only
    python average_spectrum.py output/spectra.pkl AP9MEAN.txt

    # With per-isotope line identification, 3-year mission
    python average_spectrum.py output/spectra.pkl AP9MEAN.txt --duration 3y \\
        --activities output/activities.pkl \\
        --spectra-dir result_spectra/ \\
        --save-plot avg_spectrum.pdf

    # spectra_<mode>.npz source; count rate in 20-100 keV (default: the band
    # recorded in spectra.pkl, i.e. the one given to accumulate_spectra.py)
    python average_spectrum.py output/spectra.pkl AP9MEAN.txt \\
        --activities output/activities.pkl \\
        --spectra-file result_postact/spectra_compton.npz --emin 20 --emax 100

The energy axis is that of spectra.pkl (attrs 'edges_keV'; the S-mode axis for
files written before the attr existed).  --emin/--emax only set the band of the
quoted count rate and the x limits of the plot; the spectrum keeps all channels.

Usage (library)
---------------
    from average_spectrum import compute_average_spectrum, plot_average_spectrum

    en_s, avg_spec, iso_contribs = compute_average_spectrum(
        "output/spectra.pkl", "AP9MEAN.txt",
        activities_pkl="output/activities.pkl", spectra_dir="result_spectra/",
    )
    plot_average_spectrum(en_s, avg_spec, iso_contribs, save_path="avg.pdf")
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

try:
    from activation_history import (_parse_duration, _prepare_flux_and_norm, profile_energy,
                                    belt_mask, cumulative_integral, lag_weights,
                                    step_kernel)
except ImportError:
    sys.exit("ERROR: activation_history.py not found.")

try:
    from accumulate_spectra import (load_outputs, EN_S, DELTA_E, open_spectra,
                                    band_mask, band_integral, band_label)
    from compute_activities import load_outputs as load_activities, times_of
except ImportError:
    sys.exit("ERROR: accumulate_spectra.py or compute_activities.py not found.")


# ---------------------------------------------------------------------------
# Lag-weighted average
# ---------------------------------------------------------------------------

def weighted_average(
    tau:      np.ndarray,
    values:   np.ndarray,
    dt:       float,
    w:        np.ndarray,
    w_far:    float,
    t_max:    float,
    chunk:    int = 128,
) -> np.ndarray:
    """
    sum_k w_k K_k for every column of values (n_tau x n_col), with K_k the
    step-averaged kernel: exact lags k < len(w), then the period-mean weight
    w_far for the lags up to t_max.  Returns shape (n_col,).
    """
    v = np.asarray(values, dtype=float).reshape(len(tau), -1)
    L = len(w)
    out = np.empty(v.shape[1])
    for c in range(0, v.shape[1], chunk):
        vc = v[:, c:c + chunk]
        near = w @ step_kernel(tau, vc, dt, L)
        I = cumulative_integral(tau, vc, [L * dt, max(t_max, L * dt)])
        out[c:c + chunk] = near + w_far * (I[1] - I[0]) / dt
    return out


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_average_spectrum(
    spectra_pkl:     str | Path,
    spenvis_file:    str | Path,
    activities_pkl:  str | Path | None = None,
    spectra_dir:     str | Path | None = None,
    duration_s:      float | None = None,
    all_orbit:       bool  = False,
    belt_threshold:  float = 0.0,
    spectra_file:    str | Path | None = None,
    emin:            float | None = None,
    emax:            float | None = None,
    return_info:     bool  = False,
):
    """
    Compute the steady-state out-of-belt (or, with all_orbit, orbit-averaged)
    background spectrum.

    Parameters
    ----------
    spectra_pkl     : spectra.pkl from accumulate_spectra.py
    spenvis_file    : the SPENVIS file used by accumulate_spectra.py
    activities_pkl  : activities.pkl from compute_activities.py (optional;
                      required for per-isotope decomposition)
    spectra_dir     : directory with {volume}_{isotope}_S-mode.dat files
                      (optional; with activities_pkl, for the per-isotope
                      decomposition)
    spectra_file    : or a spectra_<mode>.npz from build_spectra.py (the same
                      one as in accumulate_spectra.py)
    duration_s      : mission duration T [s] (default: end of the time grid)
    all_orbit       : average over all steps instead of the out-of-belt ones
    belt_threshold  : flux above which a step is in the belt [p/cm²/s]
    emin, emax      : band [keV] of the quoted count rate (default: the band
                      recorded in spectra.pkl, normally the whole axis)

    Returns
    -------
    en_s        : detector energy axis [keV]
    avg_spec    : average spectrum [counts/s/keV]
    iso_contribs: dict {isotope_name: avg_spectrum_array} or None
    With return_info a fourth item: dict with 'widths' [keV], 'emin', 'emax',
    'rate' (counts/s in the band), 'source', 'mode'.
    """
    spectra_pkl = Path(spectra_pkl)

    # ---- Load total spectra ------------------------------------------------
    spectra_df, _ = load_outputs(spectra_pkl.parent, fmt=spectra_pkl.suffix.lstrip("."))
    attrs = spectra_df.attrs
    if "timestep_s" not in attrs or "nprim_per_step" not in attrs:
        sys.exit(f"ERROR: {spectra_pkl} has no normalisation attrs; "
                 "rerun accumulate_spectra.py.")
    if attrs.get("activity_threshold", 0) > 0:
        print("  [WARN] spectra accumulated with an activity threshold: "
              "the average is underestimated.")

    tau   = spectra_df.index.to_numpy(dtype=float)
    S_mat = spectra_df.to_numpy(dtype=float)
    print(f"Spectra: {S_mat.shape[0]} time points × {S_mat.shape[1]} channels, "
          f"decay times {tau[0]:.2e} – {tau[-1]:.2e} s")

    # ---- Energy axis and band ------------------------------------------------
    if attrs.get("spectra_kind", "dat") == "dat" and S_mat.shape[1] == len(EN_S):
        en_s, widths = EN_S, np.full(len(EN_S), DELTA_E)
    else:
        edges = np.asarray(attrs["edges_keV"], dtype=float)
        en_s, widths = 0.5 * (edges[:-1] + edges[1:]), np.diff(edges)
    if len(en_s) != S_mat.shape[1]:
        sys.exit(f"ERROR: {spectra_pkl} has {S_mat.shape[1]} channels but its axis "
                 f"has {len(en_s)}.")
    emin = attrs.get("emin_keV") if emin is None else emin
    emax = attrs.get("emax_keV") if emax is None else emax
    mask = band_mask(en_s, emin, emax)
    src_name, src_mode = attrs.get("spectra_source", ""), attrs.get("spectra_mode", "")
    if src_name:
        print(f"Spectra source: {src_name}" + (f" (mode {src_mode})" if src_mode else ""))

    # ---- Orbit flux and lag weights ------------------------------------------
    emin = profile_energy(attrs.get("profile_emin_MeV"), None, str(spectra_pkl))
    flux, _, F_mean, dt, total = _prepare_flux_and_norm(
        spenvis_file, bool(attrs.get("in_belt_only", False)), emin)
    if abs(dt - attrs["timestep_s"]) > 1e-6 * dt:
        sys.exit(f"ERROR: time step of {spenvis_file} ({dt:.3f} s) differs from the one "
                 f"of {spectra_pkl.name} ({attrs['timestep_s']:.3f} s).")
    if duration_s is not None and duration_s > tau[-1]:
        print(f"  [WARN] duration {duration_s:.3g} s is beyond the time grid; "
              f"integrating up to {tau[-1]:.3g} s.")
    T = tau[-1] if duration_s is None else min(duration_s, tau[-1])

    select = np.ones(len(flux), bool) if all_orbit else ~belt_mask(total, belt_threshold)
    if not select.any():
        sys.exit("ERROR: no out-of-belt steps in the SPENVIS file.")
    w_far = float(np.mean(flux) / F_mean)
    w = lag_weights(flux, F_mean, select, max(1, min(len(flux), int(T // dt))))
    print(f"Mission duration T = {T:.3g} s = {T / (365.25 * 86400):.3g} yr; average over "
          + ("all steps" if all_orbit else f"the out-of-belt steps ({select.mean():.1%})"))

    def average(values: np.ndarray) -> np.ndarray:
        return weighted_average(tau, values, dt, w, w_far, T)

    # ---- Total average spectrum --------------------------------------------
    avg_spec = average(S_mat)
    rate = float(band_integral(avg_spec, widths, mask))
    print(f"Average total count rate"
          + (f" ({band_label(emin, emax)})" if band_label(emin, emax) else "")
          + f": {rate:.4e} counts/s")

    # ---- Per-isotope decomposition (optional) ------------------------------
    iso_contribs: dict[str, np.ndarray] | None = None

    if activities_pkl is not None and (spectra_dir is not None or spectra_file is not None):
        activities_pkl = Path(activities_pkl)
        src = open_spectra(spectra_dir if spectra_file is None else None, spectra_file)
        if len(src.energy) != len(en_s) or not np.allclose(src.energy, en_s, atol=1e-3):
            sys.exit(f"ERROR: the energy axis of {src.source} differs from that of "
                     f"{spectra_pkl.name}.")
        print("Computing per-isotope contributions...")

        acts_df, _ = load_activities(activities_pkl.parent,
                                     fmt=activities_pkl.suffix.lstrip("."))
        time_cols = [c for c in acts_df.columns if c.startswith("t_")]
        nprim     = {float(e): n for e, n in attrs["nprim_per_step"].items()}

        # Activity of each (volume, isotope) after one step [Bq]
        e_level = acts_df.index.get_level_values("energy_MeV").astype(float)
        keep    = e_level.isin(list(nprim))
        weights = np.array([nprim[e] for e in e_level[keep]])
        pair_act = (pd.DataFrame(acts_df.loc[keep, time_cols].to_numpy(dtype=float)
                                 * weights[:, None], index=acts_df.index[keep])
                      .groupby(level=["volume", "isotope"]).sum())
        pair_avg = average(pair_act.to_numpy().T)

        iso_acc: dict[str, np.ndarray] = {}
        for (volume, isotope), a in zip(pair_act.index, pair_avg):
            if a <= 0:
                continue
            spec = src.get(volume, isotope)
            if spec is None:
                continue
            if isotope in iso_acc:
                iso_acc[isotope] += spec * a
            else:
                iso_acc[isotope] = spec * a

        iso_contribs = iso_acc
        print(f"  {len(iso_contribs)} isotopes with spectral contributions.")

    if return_info:
        return en_s, avg_spec, iso_contribs, {
            "widths": widths, "emin": emin, "emax": emax, "rate": rate,
            "source": src_name, "mode": src_mode}
    return en_s, avg_spec, iso_contribs


# ---------------------------------------------------------------------------
# Peak finding and isotope identification
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Gamma-line database query (IAEA LiveChart REST API)
# ---------------------------------------------------------------------------

def _query_gamma_lines(
    energy_keV:  float,
    window_keV:  float = 5.0,
    cache_file:  str | Path = ".gamma_line_cache.pkl",
) -> list[str]:
    """
    Query the IAEA LiveChart REST API for nuclides that emit a gamma line
    within window_keV of energy_keV.

    Results are cached in cache_file so repeated runs do not re-query.
    Returns a list of isotope base names, e.g. ['Na22', 'I125', 'Te121m'].

    Requires internet access.  Returns [] silently on network failure.
    """
    import urllib.request, json, pickle as pkl

    cache_file = Path(cache_file)
    cache: dict = {}
    if cache_file.exists():
        try:
            with open(cache_file, "rb") as f:
                cache = pkl.load(f)
        except Exception:
            cache = {}

    key = (round(energy_keV, 1), round(window_keV, 1))
    if key in cache:
        return cache[key]

    e_lo = energy_keV - window_keV
    e_hi = energy_keV + window_keV
    # IAEA LiveChart v1 endpoint
    url = (
        "https://nds.iaea.org/relnsd/v1/data?"
        "fields=gammas&nuclides=all"
        f"&gA_lo={e_lo:.2f}&gA_hi={e_hi:.2f}"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "average_spectrum.py"})
        with urllib.request.urlopen(req, timeout=10) as r:
            raw = r.read().decode()
        # Response is CSV: z,a,symbol,...,energy,...
        lines = [l for l in raw.splitlines() if l.strip() and not l.startswith("z")]
        nuclides = set()
        for line in lines:
            parts = line.split(",")
            if len(parts) >= 3:
                sym = parts[2].strip().capitalize()
                A   = parts[1].strip()
                if sym and A.isdigit():
                    nuclides.add(f"{sym}{A}")
        result = sorted(nuclides)
    except Exception as e:
        print(f"  [WARN] gamma-line DB query failed for {energy_keV:.0f} keV: {e}")
        result = []

    cache[key] = result
    try:
        with open(cache_file, "wb") as f:
            pkl.dump(cache, f)
    except Exception:
        pass
    return result


def find_prominent_peaks(
    avg_spec:        np.ndarray,
    iso_contribs:    dict[str, np.ndarray] | None,
    n_peaks:         int   = 6,
    prominence:      float = 0.5,
    min_sep_dex:     float = 0.10,
    window_keV:      float = 5.0,
    cache_file:      str | Path = ".gamma_line_cache.pkl",
    en_s:            np.ndarray | None = None,
    band:            tuple[float | None, float | None] = (None, None),
) -> list[dict]:
    """
    Find the N tallest, well-separated peaks and identify the responsible isotope.

    Identification pipeline for each peak at energy E:
      1. Query IAEA LiveChart for nuclides with a gamma line within window_keV of E.
      2. Intersect with the isotopes present in iso_contribs (i.e. those actually
         produced in the simulation with a known spectrum file).
      3. Among the candidates, pick the one with the largest iso_contribs value
         at channel E — the one contributing most to the total spectrum.
      4. If no database match is found in iso_contribs, fall back to the largest
         iso_contribs contributor and mark it with a '?' suffix.

    Parameters
    ----------
    n_peaks      : maximum number of peaks to return
    prominence   : minimum log10-prominence noise floor
    min_sep_dex  : minimum peak separation in log10(energy) decades
    window_keV   : ±energy window for the gamma-line database query [keV]
    cache_file   : path for caching database responses between runs
    en_s         : energy axis of avg_spec [keV] (default: the S-mode axis)
    band         : (emin, emax) [keV]; only peaks inside it are considered

    Returns a list of dicts with keys:
        channel  : index into the energy axis
        energy   : keV
        value    : counts/s/keV at the peak
        isotope  : identified isotope name, 'annihilation', or None
        fraction : fraction of peak value from the identified isotope
    """
    # ---- Step 1: find peaks in total spectrum ----------------------------
    EN_S_ = EN_S if en_s is None else np.asarray(en_s)
    log_s = np.log10(np.where(avg_spec > 0, avg_spec, 1e-40))
    peaks, _ = find_peaks(log_s, prominence=prominence, width=1, distance=3)
    if len(peaks) == 0:
        return []

    if band != (None, None):
        peaks = peaks[band_mask(EN_S_, *band)[peaks]]
        if len(peaks) == 0:
            return []
    by_value = peaks[np.argsort(avg_spec[peaks])[::-1]]

    accepted: list[int] = []
    for ch in by_value:
        e = EN_S_[ch]
        too_close = any(abs(np.log10(e) - np.log10(EN_S_[a])) < min_sep_dex
                        for a in accepted)
        if not too_close:
            accepted.append(ch)
        if len(accepted) == n_peaks:
            break

    # ---- Step 2-4: identify isotope for each peak -----------------------
    results = []
    for ch in sorted(accepted):
        e_peak = EN_S_[ch]
        info: dict = {
            "channel":  ch,
            "energy":   e_peak,
            "value":    avg_spec[ch],
            "isotope":  None,
            "fraction": None,
        }

        if abs(e_peak - 511.0) < 3.0:
            info["isotope"] = "annihilation"
            results.append(info)
            continue

        if not iso_contribs:
            results.append(info)
            continue

        # Step 2: query database for candidate nuclides at this energy
        db_nuclides = _query_gamma_lines(e_peak, window_keV, cache_file)

        # Step 3: intersect with simulated isotopes.
        # iso_contribs keys may include excited states (e.g. 'Cs134[138.744]');
        # strip to base name for matching against the database.
        def base(name):
            return name.split("[")[0]

        if db_nuclides:
            candidates = {
                iso: iso_contribs[iso][ch]
                for iso in iso_contribs
                if base(iso) in db_nuclides
            }
        else:
            candidates = {}

        # Step 4: pick the largest contributor among candidates
        if candidates:
            best_iso = max(candidates, key=candidates.get)
            best_val = candidates[best_iso]
            suffix   = ""
        else:
            # No database match found among simulated isotopes
            best_iso = max(iso_contribs, key=lambda k: iso_contribs[k][ch])
            best_val = iso_contribs[best_iso][ch]
            suffix   = "?"      # uncertain attribution

        info["isotope"]  = best_iso + suffix
        info["fraction"] = best_val / avg_spec[ch] if avg_spec[ch] > 0 else 0.
        results.append(info)

    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_spectrum(en_s, avg_spec, path):
    np.savetxt(path,
               np.column_stack([en_s, avg_spec]),
               header="energy_keV   avg_rate_counts_per_s_per_keV",
               fmt="%.4f  %.6e")
    print(f"Average spectrum saved to {path}")


# ---------------------------------------------------------------------------
# Publication-quality plot
# ---------------------------------------------------------------------------

def plot_average_spectrum(
    en_s:         np.ndarray,
    avg_spec:     np.ndarray,
    iso_contribs: dict[str, np.ndarray] | None = None,
    n_label:      int   = 5,
    prominence:   float = 0.5,
    min_sep_dex:  float = 0.10,
    window_keV:   float = 5.0,
    cache_file:   str | Path = ".gamma_line_cache.pkl",
    title:        str   = "Steady-state out-of-belt background spectrum",
    save_path:    str | Path | None = None,
    widths:       np.ndarray | None = None,
    emin:         float | None = None,
    emax:         float | None = None,
    xscale:       str   = "log",
) -> None:
    """
    Publication-quality log-log plot of the orbit-averaged spectrum.

    The N most prominent spectral peaks are labelled with the dominant
    isotope name and energy in keV.  An arrow connects the label to the peak.

    Figure is sized for a single journal column (3.5 in wide).

    widths     : channel widths [keV] for the quoted count rate (default: the
                 2 keV of the S-mode axis)
    emin, emax : energy band [keV] of the count rate, also the x limits and the
                 range where peaks are labelled (default: the whole axis)
    xscale     : 'log' or 'linear' energy axis
    """
    peaks_info = find_prominent_peaks(
        avg_spec, iso_contribs, n_peaks=n_label,
        prominence=prominence, min_sep_dex=min_sep_dex,
        window_keV=window_keV, cache_file=cache_file,
        en_s=en_s, band=(emin, emax),
    )

    pos        = avg_spec > 0
    w          = np.full(len(en_s), DELTA_E) if widths is None else np.asarray(widths)
    total_rate = float(band_integral(avg_spec, w, band_mask(en_s, emin, emax)))
    blab       = band_label(emin, emax)

    # ---- Figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    # Spectrum
    ax.plot(en_s[pos], avg_spec[pos],
            drawstyle="steps-mid",
            color="#2166ac", linewidth=0.9, zorder=3,
            label=rf"$\langle S \rangle$  ({total_rate:.2e} cts/s"
                  + (f", {blab}" if blab else "") + ")")

    # Per-isotope contributions as lighter fills (optional, if available)
    if iso_contribs and len(peaks_info) > 0:
        # Identify the unique dominant isotopes across labelled peaks
        dominant = list({p["isotope"] for p in peaks_info if p["isotope"]})
        cmap     = plt.colormaps["tab10"]
        iso_colors: dict[str, str] = {
            iso: cmap(i % 10) for i, iso in enumerate(dominant)
        }
        for iso in dominant:
            if iso not in iso_contribs:   # e.g. 'annihilation' pseudo-label
                continue
            s = iso_contribs[iso]
            pos_i = s > 0
            if pos_i.any():
                ax.fill_between(en_s[pos_i], s[pos_i], 0,
                                step="mid", alpha=0.18,
                                color=iso_colors[iso], zorder=2)

    # ---- Peak annotations --------------------------------------------------
    y_lo, y_hi = avg_spec[pos].min(), avg_spec[pos].max()
    log_range  = np.log10(y_hi) - np.log10(y_lo)

    # Sort peaks by energy for left-to-right annotation.
    # Alternate label positions above/below the spectrum to reduce overlap.
    peaks_info_sorted = sorted(peaks_info, key=lambda p: p["energy"])

    ax.set_xscale(xscale)
    ax.set_yscale("log")
    if emin is not None or emax is not None:
        ax.set_xlim(emin if emin is not None else (max(en_s[0], 1.0) if xscale == "log"
                                                   else en_s[0]),
                    emax if emax is not None else en_s[-1])

    # Two rows of label positions in log-y space: upper and lower band
    upper_fracs = np.linspace(0.80, 0.95, (len(peaks_info_sorted) + 1) // 2)
    lower_fracs = np.linspace(0.60, 0.75, len(peaks_info_sorted) // 2)
    # interleave: even indices -> upper, odd -> lower
    y_fracs = []
    ui = li = 0
    for k in range(len(peaks_info_sorted)):
        if k % 2 == 0 and ui < len(upper_fracs):
            y_fracs.append(upper_fracs[ui]); ui += 1
        elif li < len(lower_fracs):
            y_fracs.append(lower_fracs[li]); li += 1
        else:
            y_fracs.append(upper_fracs[min(ui, len(upper_fracs)-1)]); ui += 1

    for k, (pk, y_frac) in enumerate(zip(peaks_info_sorted, y_fracs)):
        e_peak = pk["energy"]
        v_peak = pk["value"]
        iso    = pk["isotope"]

        # Label: isotope name on first line, energy on second
        if iso:
            label = f"{iso}\n{e_peak:.0f} keV"
        else:
            label = f"{e_peak:.0f} keV"

        y_text = 10 ** (np.log10(y_lo) + y_frac * log_range)

        ax.annotate(
            label,
            xy=(e_peak, v_peak),
            xytext=(e_peak, y_text),
            fontsize=5,
            ha="center", va="bottom",
            color="#333333",
            arrowprops=dict(
                arrowstyle="-",
                color="#999999",
                linewidth=0.6,
                shrinkA=2, shrinkB=2,
            ),
            annotation_clip=False,
        )

    # ---- Axes formatting ---------------------------------------------------
    ax.set_xlabel("Energy (keV)", fontsize=8)
    ax.set_ylabel(r"Rate (cts s$^{-1}$ keV$^{-1}$)", fontsize=8)
    ax.tick_params(axis="both", which="both",
                   direction="in", top=True, right=True, labelsize=7)
    ax.grid(True, which="major", linestyle=":", linewidth=0.4,
            color="gray", alpha=0.5)
    ax.set_title(title, fontsize=8, pad=4)
    ax.legend(fontsize=6.5, framealpha=0.8, loc="lower left")

    fig.tight_layout(pad=0.5)

    if save_path is not None:
        fig.savefig(Path(save_path), dpi=300, bbox_inches="tight")
        print(f"Spectrum plot saved to {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Compute the steady-state orbit-averaged background spectrum.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("spectra_pkl",     help="spectra.pkl from accumulate_spectra.py")
    p.add_argument("spenvis_file",    help="SPENVIS file used by accumulate_spectra.py")
    p.add_argument("--duration",      default=None, metavar="T",
                   help="Mission duration, e.g. 3y, 18mo (default: end of the time grid).")
    p.add_argument("--all-orbit",     action="store_true",
                   help="Average over all time steps instead of the out-of-belt ones.")
    p.add_argument("--belt-threshold", type=float, default=0.0, metavar="FLUX",
                   help="Total flux [p/cm²/s] above which a step is in the belt.")
    p.add_argument("--activities",    default=None, metavar="PKL",
                   help="activities.pkl for per-isotope line identification.")
    p.add_argument("--spectra-dir",   default=None, metavar="DIR",
                   help="result_spectra/ directory (.dat files) for the per-isotope "
                        "decomposition.")
    p.add_argument("--spectra-file",  default=None, metavar="NPZ",
                   help="spectra_<mode>.npz from build_spectra.py instead of "
                        "--spectra-dir (the one used in accumulate_spectra.py).")
    p.add_argument("--emin", type=float, default=None, metavar="KEV",
                   help="Lower edge of the energy band of the count rate and plot "
                        "[keV] (default: the band recorded in spectra.pkl).")
    p.add_argument("--emax", type=float, default=None, metavar="KEV",
                   help="Upper edge of the energy band [keV].")
    p.add_argument("--xscale", choices=["log", "linear"], default="log",
                   help="Energy axis of the plot.")
    p.add_argument("--n-label",   type=int, default=5, metavar="N",
                   help="Number of peaks to label.")
    p.add_argument("--prominence", type=float, default=0.5,
                   help="Minimum log10-prominence for peak detection.")
    p.add_argument("--min-sep", type=float, default=0.10, metavar="DEX",
                   help="Minimum separation between labelled peaks in decades (log10 of energy ratio). Default 0.10 = factor ~1.26 in energy.")
    p.add_argument("--window", type=float, default=5.0, metavar="KEV",
                   help="+/- energy window for gamma-line DB query [keV].")
    p.add_argument("--cache", default=".gamma_line_cache.pkl", metavar="FILE",
                   help="Cache file for gamma-line DB responses.")
    p.add_argument("--save-dat",  default=None, metavar="PATH")
    p.add_argument("--save-plot", default=None, metavar="PATH")
    args = p.parse_args(argv)

    en_s, avg_spec, iso_contribs, info = compute_average_spectrum(
        spectra_pkl=args.spectra_pkl,
        spenvis_file=args.spenvis_file,
        activities_pkl=args.activities,
        spectra_dir=args.spectra_dir,
        duration_s=_parse_duration(args.duration) if args.duration else None,
        all_orbit=args.all_orbit,
        belt_threshold=args.belt_threshold,
        spectra_file=args.spectra_file,
        emin=args.emin,
        emax=args.emax,
        return_info=True,
    )

    if args.save_dat:
        save_spectrum(en_s, avg_spec, args.save_dat)

    plot_average_spectrum(
        en_s, avg_spec, iso_contribs,
        n_label=args.n_label,
        prominence=args.prominence,
        min_sep_dex=args.min_sep,
        window_keV=args.window,
        cache_file=args.cache,
        save_path=args.save_plot,
        widths=info["widths"],
        emin=info["emin"],
        emax=info["emax"],
        xscale=args.xscale,
    )


if __name__ == "__main__":
    main()
