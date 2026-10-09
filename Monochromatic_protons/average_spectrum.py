"""
average_spectrum.py
===================
Compute the steady-state orbit-averaged background spectrum after a long
mission duration, and identify the isotopes responsible for the most
prominent spectral lines.

Physical basis
--------------
After many orbits the activation has reached a periodic steady state.
The orbit-averaged spectrum is:

    <S(E)> = (1 / T_orb) * integral_0^inf S(tau, E) dtau

where S(tau, E) [counts/s/keV] is the spectrum at decay time tau after
a single reference irradiation (from accumulate_spectra.py).

Per-isotope decomposition
-------------------------
If activities_pkl and spectra_dir are supplied, the contribution of each
isotope i (summed over all volumes and simulation energies) is:

    <S_i(E)> = (1 / T_orb) * sum_{v,e} spec_{i,v}(E)
               * n_prim(e) * integral_0^inf A_{i,v,e}(tau) dtau

where spec_{i,v} is the background spectrum file [counts/s/keV/decay],
n_prim(e) is the number of primaries for energy band e, and A_{i,v,e}(tau)
is the normalised activity [Bq/primary].  The dominant isotope at each
spectral peak is then the one with the largest <S_i(E_peak)>.

Usage (CLI)
-----------
    # Total spectrum only
    python average_spectrum.py spectra.pkl count_rate.dat AP9MEAN.txt

    # With per-isotope line identification
    python average_spectrum.py spectra.pkl count_rate.dat AP9MEAN.txt \\
        --activities activities.pkl \\
        --spectra-dir result_spectra/ \\
        --save-plot avg_spectrum.pdf

Usage (library)
---------------
    from average_spectrum import compute_average_spectrum, plot_average_spectrum

    en_s, avg_spec, iso_contribs = compute_average_spectrum(
        "spectra.pkl", "count_rate.dat", "AP9MEAN.txt",
        activities_pkl="activities.pkl", spectra_dir="result_spectra/",
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

# Compatibility: np.trapezoid (NumPy >= 2.0) vs np.trapz (NumPy < 2.0)
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)

try:
    from activation_history import _read_total_flux
except ImportError:
    sys.exit("ERROR: activation_history.py not found.")

try:
    from accumulate_spectra import load_outputs, EN_S, load_spectrum, spectrum_path
    from compute_activities import load_outputs as load_activities
except ImportError:
    sys.exit("ERROR: accumulate_spectra.py or compute_activities.py not found.")

try:
    from spenvis_parser import parse_spenvis, band_fluxes_to_primaries
except ImportError:
    sys.exit("ERROR: spenvis_parser.py not found.")


# ---------------------------------------------------------------------------
# Integration helper
# ---------------------------------------------------------------------------

def _integrate_rows(
    tau:    np.ndarray,
    values: np.ndarray,
    n_interp: int = 1000,
) -> np.ndarray:
    """
    Integrate rows of `values` (shape n_tau × n_col) over tau using a fine
    log-spaced grid.  Returns shape (n_col,).

    Uses lin-log interpolation (linear in tau, log in value) for accuracy
    on the coarse decade-spaced tau grid.
    """
    tau_fine  = np.logspace(np.log10(tau[0]), np.log10(tau[-1]), n_interp)
    n_col     = values.shape[1] if values.ndim == 2 else 1
    vals_2d   = values.reshape(len(tau), -1)
    fine      = np.zeros((n_interp, vals_2d.shape[1]))

    for c in range(vals_2d.shape[1]):
        col = vals_2d[:, c]
        pos = col > 0
        if not pos.any():
            continue
        tau_p = tau[pos]
        log_v = np.log(col[pos])
        inside = (tau_fine >= tau_p[0]) & (tau_fine <= tau_p[-1])
        if inside.any():
            fine[inside, c] = np.exp(np.interp(tau_fine[inside], tau_p, log_v))

    return _trapz(fine, tau_fine, axis=0)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_average_spectrum(
    spectra_pkl:     str | Path,
    count_rate_file: str | Path,
    spenvis_file:    str | Path,
    activities_pkl:  str | Path | None = None,
    spectra_dir:     str | Path | None = None,
    energies:        list[float] | None = None,
    timestep:        float = 60.0,
    R:               float = 5000.0,
    thetamax:        float = 0.8021,
    n_interp:        int   = 1000,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray] | None]:
    """
    Compute the steady-state orbit-averaged background spectrum.

    Parameters
    ----------
    spectra_pkl     : spectra.pkl from accumulate_spectra.py
    count_rate_file : count_rate.dat from accumulate_spectra.py
    spenvis_file    : SPENVIS AP9/AE9 file (for T_orb)
    activities_pkl  : activities.pkl from compute_activities.py (optional;
                      required for per-isotope decomposition)
    spectra_dir     : directory with {volume}_{isotope}_S-mode.dat files
                      (optional; required for per-isotope decomposition)
    energies        : band edge energies [MeV] matching simulation grid
    timestep, R, thetamax : simulation geometry (for n_prim calculation)
    n_interp        : points in the fine tau integration grid

    Returns
    -------
    en_s        : detector energy axis [keV]
    avg_spec    : orbit-averaged spectrum [counts/s/keV]
    iso_contribs: dict {isotope_name: avg_spectrum_array} or None
    """
    if energies is None:
        energies = [7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 150, 200, 300, 400]

    spectra_pkl     = Path(spectra_pkl)
    count_rate_file = Path(count_rate_file)

    # ---- Load total spectra ------------------------------------------------
    spectra_df, _ = load_outputs(
        str(spectra_pkl.parent), fmt=spectra_pkl.suffix.lstrip(".")
    )
    tau   = spectra_df.index.to_numpy(dtype=float)
    S_mat = spectra_df.values.astype(float)
    print(f"Spectra: {S_mat.shape[0]} time points × {S_mat.shape[1]} channels")
    print(f"Decay time range: {tau[0]:.2e} – {tau[-1]:.2e} s")

    # ---- T_orb from SPENVIS ------------------------------------------------
    _, mjd = _read_total_flux(spenvis_file)
    dt     = float(np.median(np.diff(mjd)) * 86400.0)
    T_orb  = (mjd[-1] - mjd[0]) * 86400.0 + dt
    print(f"Orbital period T_orb = {T_orb:.1f} s = {T_orb/3600:.2f} h")

    # ---- Total average spectrum --------------------------------------------
    integral = _integrate_rows(tau, S_mat, n_interp)
    avg_spec = integral / T_orb
    total_rate = avg_spec.sum() * 2.0   # dE = 2 keV
    print(f"Orbit-averaged total count rate: {total_rate:.4e} counts/s")

    # ---- Per-isotope decomposition (optional) ------------------------------
    iso_contribs: dict[str, np.ndarray] | None = None

    if activities_pkl is not None and spectra_dir is not None:
        activities_pkl = Path(activities_pkl)
        spectra_dir    = Path(spectra_dir)
        print("Computing per-isotope contributions...")

        acts_df, _ = load_activities(
            str(activities_pkl.parent), fmt=activities_pkl.suffix.lstrip(".")
        )
        time_cols = [c for c in acts_df.columns if c.startswith("t_")]
        act_tau   = np.array([float(c[2:]) for c in time_cols])

        # n_primaries per band (same calculation as in accumulate_spectra.py)
        _, band_fluxes = parse_spenvis(spenvis_file, energies=energies)
        n_prim_per_band = band_fluxes_to_primaries(
            energies, band_fluxes, timestep=timestep, R=R, thetamax=thetamax
        )
        sim_energies  = energies[:-1]
        e_to_nprim    = {float(e): n for e, n in zip(sim_energies, n_prim_per_band)}

        # Accumulate per-isotope: sum contributions across all volumes and energies
        iso_acc: dict[str, np.ndarray] = {}
        n_channels = len(EN_S)

        for (energy_MeV, volume, isotope), row in acts_df.iterrows():
            n_prim = e_to_nprim.get(energy_MeV)
            if n_prim is None:
                continue

            spec_file = spectrum_path(spectra_dir, volume, isotope)
            spec = load_spectrum(spec_file)
            if spec is None:
                continue

            # Integral of A(tau) dtau [decays/primary * s] on the activity grid
            A_arr = row[time_cols].to_numpy(dtype=float)
            pos   = A_arr > 0
            if not pos.any():
                continue

            # Integrate using the fine grid
            tau_p  = act_tau[pos]
            A_p    = A_arr[pos]
            if len(tau_p) < 2:
                continue
            tau_fine_i = np.logspace(np.log10(tau_p[0]), np.log10(tau_p[-1]),
                                     min(n_interp, 500))
            log_A   = np.log(A_p)
            inside  = (tau_fine_i >= tau_p[0]) & (tau_fine_i <= tau_p[-1])
            A_fine  = np.zeros(len(tau_fine_i))
            A_fine[inside] = np.exp(np.interp(tau_fine_i[inside], tau_p, log_A))
            int_A   = _trapz(A_fine, tau_fine_i)   # [Bq/primary * s = decays/primary]

            # Contribution to average spectrum [counts/s/keV]
            contrib = spec * (n_prim * int_A / T_orb)

            if isotope in iso_acc:
                iso_acc[isotope] += contrib
            else:
                iso_acc[isotope] = contrib.copy()

        iso_contribs = iso_acc
        print(f"  {len(iso_contribs)} isotopes with spectral contributions.")

    return EN_S, avg_spec, iso_contribs


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

    Returns a list of dicts with keys:
        channel  : index into EN_S
        energy   : keV
        value    : counts/s/keV at the peak
        isotope  : identified isotope name, 'annihilation', or None
        fraction : fraction of peak value from the identified isotope
    """
    # ---- Step 1: find peaks in total spectrum ----------------------------
    log_s = np.log10(np.where(avg_spec > 0, avg_spec, 1e-40))
    peaks, _ = find_peaks(log_s, prominence=prominence, width=1, distance=3)
    if len(peaks) == 0:
        return []

    by_value = peaks[np.argsort(avg_spec[peaks])[::-1]]

    accepted: list[int] = []
    for ch in by_value:
        e = EN_S[ch]
        too_close = any(abs(np.log10(e) - np.log10(EN_S[a])) < min_sep_dex
                        for a in accepted)
        if not too_close:
            accepted.append(ch)
        if len(accepted) == n_peaks:
            break

    # ---- Step 2-4: identify isotope for each peak -----------------------
    results = []
    for ch in sorted(accepted):
        e_peak = EN_S[ch]
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
    title:        str   = "Steady-state orbit-averaged background spectrum",
    save_path:    str | Path | None = None,
) -> None:
    """
    Publication-quality log-log plot of the orbit-averaged spectrum.

    The N most prominent spectral peaks are labelled with the dominant
    isotope name and energy in keV.  An arrow connects the label to the peak.

    Figure is sized for a single journal column (3.5 in wide).
    """
    peaks_info = find_prominent_peaks(
        avg_spec, iso_contribs, n_peaks=n_label,
        prominence=prominence, min_sep_dex=min_sep_dex,
        window_keV=window_keV, cache_file=cache_file,
    )

    pos        = avg_spec > 0
    total_rate = avg_spec.sum() * 2.0   # dE = 2 keV

    # ---- Figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    # Spectrum
    ax.plot(en_s[pos], avg_spec[pos],
            drawstyle="steps-mid",
            color="#2166ac", linewidth=0.9, zorder=3,
            label=rf"$\langle S \rangle$  ({total_rate:.2e} cts/s)")

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

    ax.set_xscale("log")
    ax.set_yscale("log")

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
    p.add_argument("count_rate_file", help="count_rate.dat from accumulate_spectra.py")
    p.add_argument("spenvis_file",    help="SPENVIS AP9/AE9 output file")
    p.add_argument("--activities",    default=None, metavar="PKL",
                   help="activities.pkl for per-isotope line identification.")
    p.add_argument("--spectra-dir",   default=None, metavar="DIR",
                   help="result_spectra/ directory for per-isotope decomposition.")
    p.add_argument("--energies", nargs="+", type=float, metavar="E",
                   default=[7,10,15,20,30,40,50,60,70,100,150,200,300,400])
    p.add_argument("--timestep",  type=float, default=60.0)
    p.add_argument("--R",         type=float, default=5000.0)
    p.add_argument("--thetamax",  type=float, default=0.8021)
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
    p.add_argument("--n-interp",  type=int, default=1000)
    p.add_argument("--save-dat",  default=None, metavar="PATH")
    p.add_argument("--save-plot", default=None, metavar="PATH")
    args = p.parse_args(argv)

    en_s, avg_spec, iso_contribs = compute_average_spectrum(
        spectra_pkl=args.spectra_pkl,
        count_rate_file=args.count_rate_file,
        spenvis_file=args.spenvis_file,
        activities_pkl=args.activities,
        spectra_dir=args.spectra_dir,
        energies=args.energies,
        timestep=args.timestep,
        R=args.R,
        thetamax=args.thetamax,
        n_interp=args.n_interp,
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
    )


if __name__ == "__main__":
    main()
