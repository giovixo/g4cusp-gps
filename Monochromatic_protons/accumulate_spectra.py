"""
accumulate_spectra.py
=====================
Background spectrum and count rate at each decay time after one irradiation
step at the orbit-mean proton spectrum: the reference curve R(tau) that
activation_history.py convolves with the orbital flux, and that
average_spectrum.py integrates.

Method (Campana et al. 2026, Eqs. 4 and 6)
------------------------------------------
The orbit-mean proton fluence of one SPENVIS time step dt is assigned to the
simulation energies E_j, each standing for the band [e_j, e_{j+1}] around it
(edges at the geometric midpoints between simulation energies; the outer edges
are the first and last energy).  The primaries of energy E_j in one step are

    N_j = dt * (F(>e_j) - F(>e_{j+1})) * A_beam            (s_j dE_j of Eq. 4)

where F(>E) is the orbit-mean integral flux and A_beam = pi R^2 sin^2(thetamax)
is the beam area of the simulation source (from run_info.json, carried in the
activities.pkl attrs).  The activity and the spectrum are then

    A_{v,i}(tau) = sum_j N_j a_{j,v,i}(tau)                [Bq]
    S(tau, E)    = sum_{v,i} A_{v,i}(tau) F_{v,i}(E)       [counts/s/keV]

with a the normalised activity [Bq/primary] from compute_activities.py and
F_{v,i} the detector spectrum per decay [counts/keV/decay].

Inputs
------
  activities.pkl      from compute_activities.py
  SPENVIS file        AP9/AE9 integral flux along the orbit (dt and mean flux)
  detector spectra per decay [counts/keV/decay], from one of two sources:
    --spectra-dir DIR    result_spectra/{volume}_{isotope}_S-mode.dat, one value
                         per channel of the fixed S-mode axis (2047 channels,
                         18-4112 keV, 2 keV wide) -- the default;
    --spectra-file FILE  spectra_<mode>.npz from build_spectra.py (modes
                         scat_single, abs_single, compton, any); the energy axis
                         (default 0-2000 keV, 1 keV bins) comes from the file.

Output
------
  spectra.pkl / .parquet
      DataFrame: index time_s, columns = detector channel centres [keV].
      Spectrum [counts/s/keV] at each decay time.  attrs: 'timestep_s', 'duty',
      'beam_area_cm2', 'nprim_per_step' ({E: N_j}), 'in_belt_only',
      'activity_threshold', and the spectra source: 'spectra_source' (file or
      directory), 'spectra_mode', 'edges_keV' (channel edges), 'emin_keV',
      'emax_keV' (band of the count rate; None = whole axis).
  count_rate.dat
      Two-column ASCII: time_s  count_rate [counts/s], the spectrum integrated
      over the band (sum of spectrum x channel width over the channels whose
      centres lie in [emin, emax]; default the whole axis).  With a spectra file
      or a band, header lines 'source:', 'mode:' and 'band_keV:' record them.

Usage (CLI)
-----------
    python accumulate_spectra.py output/activities.pkl AP9MEAN.txt
    python accumulate_spectra.py output/activities.pkl AP9MEAN.txt \\
        --spectra-dir result_spectra/ --outdir output/ --save-plot count_rate.pdf

    # spectra_<mode>.npz from build_spectra.py, count rate in 20-100 keV
    python accumulate_spectra.py output/activities.pkl AP9MEAN.txt \\
        --spectra-file result_postact/spectra_compton.npz --emin 20 --emax 100 \\
        --outdir output_compton/

Usage (library)
---------------
    from accumulate_spectra import accumulate

    spectra_df, count_rate = accumulate(
        "output/activities.pkl", "AP9MEAN.txt",
        spectra_dir="result_spectra/",          # or spectra_file="spectra_compton.npz"
        emin=20.0, emax=100.0,                  # optional band of the count rate
    )
"""

from __future__ import annotations

import argparse
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
    from compute_activities import load_outputs as load_activities, times_of
except ImportError:
    sys.exit("ERROR: compute_activities.py not found.")

try:
    from spenvis_parser import parse_spenvis
except ImportError:
    sys.exit("ERROR: spenvis_parser.py not found.")

try:
    from activation_history import _read_total_flux
except ImportError:
    sys.exit("ERROR: activation_history.py not found.")


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
DELTA_E = 2.0            # keV per channel


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
# Spectra source: the .dat files (fixed S-mode axis) or a spectra_<mode>.npz
# ---------------------------------------------------------------------------

class SpectraSource:
    """
    Detector spectra per decay of the (volume, isotope) pairs, from either
    source, with the energy axis they are tabulated on.

    energy  : channel centres [keV]
    widths  : channel widths [keV]
    edges   : channel edges [keV] (len(energy) + 1)
    get(volume, isotope) -> counts/keV/decay array, or None if absent
    kind    : 'dat' or 'npz';  source : path of the file or directory;
    mode    : event class of an npz file, '' for the .dat files
    """

    def __init__(self, spectra_dir=None, spectra_file=None):
        if spectra_file is not None and spectra_dir is not None:
            raise ValueError("Give either spectra_dir or spectra_file, not both.")
        if spectra_file is not None:
            from pair_spectra import load_pair_spectra
            self._lib = load_pair_spectra(spectra_file)
            self.kind, self.source = "npz", str(spectra_file)
            self.mode = self._lib.mode
            self.energy = self._lib.centres
            self.widths = self._lib.widths
            self.edges = self._lib.edges
            self._get = self._lib.spectrum
        else:
            d = Path("result_spectra" if spectra_dir is None else spectra_dir)
            self.kind, self.source, self.mode = "dat", str(d), ""
            self.energy = EN_S
            self.widths = np.full(N_CHANNELS, DELTA_E)
            self.edges = np.concatenate([EN_S - DELTA_E / 2, [EN_S[-1] + DELTA_E / 2]])
            self._get = lambda v, i: load_spectrum(spectrum_path(d, v, i))

    def get(self, volume: str, isotope: str) -> np.ndarray | None:
        return self._get(volume, isotope)

    @property
    def default(self) -> bool:
        """True for the .dat source (the original behaviour)."""
        return self.kind == "dat"


def open_spectra(spectra_dir=None, spectra_file=None) -> SpectraSource:
    """Spectra source from --spectra-dir (.dat files) or --spectra-file (.npz)."""
    return SpectraSource(spectra_dir, spectra_file)


def band_mask(energy: np.ndarray, emin: float | None = None,
              emax: float | None = None) -> np.ndarray:
    """Channels whose centres lie in [emin, emax] keV (None = no limit)."""
    m = np.ones(len(energy), bool)
    if emin is not None:
        m &= energy >= emin
    if emax is not None:
        m &= energy <= emax
    if not m.any():
        raise ValueError(f"No channel centre in the band [{emin}, {emax}] keV "
                         f"(axis {energy[0]:g}-{energy[-1]:g} keV).")
    return m


def band_integral(spec: np.ndarray, widths: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """
    Integral over the channels in mask of spec (last axis, counts/s/keV) times
    the channel widths: counts/s.  Uniform widths and the whole axis give the
    plain sum times the width.
    """
    spec = np.asarray(spec)
    if mask.all() and np.all(widths == widths[0]):
        return spec.sum(axis=-1) * widths[0]
    return spec[..., mask] @ widths[mask]


def band_label(emin: float | None, emax: float | None) -> str:
    """'20-100 keV', '>= 20 keV', ... or '' for the whole axis."""
    if emin is None and emax is None:
        return ""
    if emin is not None and emax is not None:
        return f"{emin:g}-{emax:g} keV"
    return f">= {emin:g} keV" if emin is not None else f"<= {emax:g} keV"


# ---------------------------------------------------------------------------
# Primaries per orbital time step
# ---------------------------------------------------------------------------

def band_edges(energies: list[float] | np.ndarray) -> np.ndarray:
    """
    Edges of the bands represented by the simulation energies: geometric
    midpoints between consecutive energies, and the first and last energy as
    outer edges.  Band j = [edges[j], edges[j+1]] contains energies[j].
    """
    e = np.sort(np.asarray(energies, dtype=float))
    if len(e) < 2:
        raise ValueError("At least two simulation energies are needed to define the bands.")
    return np.concatenate([[e[0]], np.sqrt(e[:-1] * e[1:]), [e[-1]]])


def source_beam_area(
    attrs:    dict,
    R:        float | None = None,
    thetamax: float | None = None,
) -> float:
    """
    Beam area pi R^2 sin^2(thetamax) [cm^2] of the simulation source: from R and
    thetamax if given, otherwise from the source geometry in the attrs of
    activities.pkl (run_info.json of 0_run.py).
    """
    if R is not None or thetamax is not None:
        if R is None or thetamax is None:
            raise ValueError("Give both R and thetamax, or neither.")
        return float(np.pi * R**2 * np.sin(np.radians(thetamax))**2)
    source = attrs.get("source") or {}
    if "beam_area_cm2" not in source:
        raise ValueError("The activities carry no source geometry (run_info.json "
                         "missing at step 1): give R and thetamax.")
    return float(source["beam_area_cm2"])


def tabulated_energies(spenvis_file: str | Path) -> np.ndarray:
    """Energy levels [MeV] of a SPENVIS integral flux file (header line)."""
    with open(spenvis_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if "Energy levels" in line:
                return np.array([float(x) for x in line.split(":", 1)[1].split()])
    raise ValueError(f"No 'Energy levels' line in {spenvis_file}")


def primaries_per_step(
    spenvis_file:  str | Path,
    energies:      list[float],
    beam_area_cm2: float,
    in_belt_only:  bool = False,
) -> tuple[dict[float, float], float, float]:
    """
    Primaries of each simulation energy in one SPENVIS time step at the mean
    proton spectrum (see the module docstring).

    Returns
    -------
    nprim    : {energy_MeV: primaries per time step}
    timestep : SPENVIS time step dt [s]
    duty     : orbit mean of F(t) / F_mean, the factor that turns the response
               to one step at F_mean into the orbit average: 1, or the in-belt
               fraction of the time steps if in_belt_only
    """
    flux, mjd = _read_total_flux(spenvis_file)
    timestep = float(np.median(np.diff(mjd)) * 86400.0)
    duty = float((flux > 0).mean()) if in_belt_only else 1.0

    energies = sorted(float(e) for e in energies)
    edges = band_edges(energies)
    e_max = tabulated_energies(spenvis_file)[-1]
    if edges[-1] > e_max:
        warnings.warn(f"{Path(spenvis_file).name} tabulates the flux up to {e_max:g} MeV: "
                      f"band edges above it are clipped, so protons above {e_max:g} MeV are "
                      f"ignored and the simulation energies above it get no primaries.",
                      stacklevel=2)
        edges = np.minimum(edges, e_max)
    _, band_fluxes = parse_spenvis(spenvis_file, energies=list(edges),
                                   in_belt_only=in_belt_only)
    nprim = timestep * np.asarray(band_fluxes) * beam_area_cm2
    return dict(zip(energies, nprim.tolist())), timestep, duty


# ---------------------------------------------------------------------------
# Core accumulation
# ---------------------------------------------------------------------------

def accumulate(
    activities_pkl:     str | Path,
    spenvis_file:       str | Path,
    spectra_dir:        str | Path  = "result_spectra",
    energies:           list[float] | None = None,
    activity_threshold: float       = 0.0,
    R:                  float | None = None,
    thetamax:           float | None = None,
    in_belt_only:       bool        = False,
    spectra_file:       str | Path | None = None,
    emin:               float | None = None,
    emax:               float | None = None,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Background spectrum and count rate at each decay time after one SPENVIS
    time step of irradiation at the orbit-mean proton spectrum.

    Parameters
    ----------
    activities_pkl     : activities.pkl from compute_activities.py
    spenvis_file       : SPENVIS AP9/AE9 integral flux file
    spectra_dir        : directory with the {volume}_{isotope}_S-mode.dat files
    spectra_file       : alternatively a spectra_<mode>.npz from build_spectra.py
                         (then spectra_dir is ignored and the energy axis, the
                         channel widths and the spectra come from the file)
    emin, emax         : band [keV] of the count rate: the spectrum integrated
                         over the channels whose centres lie in it (default: the
                         whole axis).  The spectra keep all channels.
    energies           : simulation energies to use [MeV] (default: all the
                         energies of activities.pkl)
    activity_threshold : (volume, isotope) pairs are skipped at the decay times
                         where their activity after the single step is not
                         above this value [Bq].  The default 0 keeps everything;
                         a cut on the single-step activity drops contributions
                         that add up over many steps in activation_history.py
                         and average_spectrum.py, so use it only as a speed-up.
    R, thetamax        : source sphere radius [cm] and cone half-angle [deg];
                         by default the beam area comes from the activities
    in_belt_only       : mean flux over the in-belt time steps only (use the
                         same choice in activation_history.py)

    Returns
    -------
    spectra_df : DataFrame, index time_s, columns detector channels [keV];
                 spectrum [counts/s/keV].  The attrs record the normalisation
                 (see the module docstring).
    count_rate : Series, index time_s; total count rate [counts/s].
    """
    activities_pkl = Path(activities_pkl)
    src  = open_spectra(None if spectra_file is not None else spectra_dir, spectra_file)
    mask = band_mask(src.energy, emin, emax)
    print(f"Spectra: {src.source}" + (f" (mode {src.mode})" if src.mode else "")
          + f", {len(src.energy)} channels {src.edges[0]:g}-{src.edges[-1]:g} keV"
          + (f"; count rate band {band_label(emin, emax)}" if band_label(emin, emax) else ""))

    # ---- Normalised activities ----------------------------------------
    print(f"Loading activities from: {activities_pkl}")
    activities_df, _ = load_activities(activities_pkl.parent,
                                       fmt=activities_pkl.suffix.lstrip("."))
    times_s   = times_of(activities_df)
    time_cols = [c for c in activities_df.columns if c.startswith("t_")]
    print(f"  {len(activities_df)} (energy, volume, isotope) entries, "
          f"{len(times_s)} times: {times_s[0]:.2e} – {times_s[-1]:.2e} s")

    sim_energies = sorted(float(e) for e in
                          activities_df.index.get_level_values("energy_MeV").unique())
    if energies is not None:
        unknown = sorted(set(map(float, energies)) - set(sim_energies))
        if unknown:
            raise ValueError(f"Energies {unknown} MeV are not in {activities_pkl.name} "
                             f"(simulated: {sim_energies}).")
        sim_energies = sorted(map(float, energies))

    # ---- Primaries per time step ---------------------------------------
    beam_area = source_beam_area(activities_df.attrs, R, thetamax)
    print(f"Loading SPENVIS flux from: {spenvis_file}")
    nprim, timestep, duty = primaries_per_step(spenvis_file, sim_energies, beam_area,
                                               in_belt_only=in_belt_only)
    edges = band_edges(sim_energies)
    print(f"  time step {timestep:.1f} s, beam area {beam_area:.4g} cm^2"
          + (f", in-belt fraction {duty:.3f}" if in_belt_only else ""))
    for j, e in enumerate(sim_energies):
        print(f"    {e:>7.1f} MeV  [{edges[j]:7.2f}, {edges[j+1]:7.2f}]: "
              f"{nprim[e]:.4e} primaries/step")

    # ---- Activity of each (volume, isotope) after one step [Bq] --------
    e_level = activities_df.index.get_level_values("energy_MeV").astype(float)
    keep    = e_level.isin(sim_energies)
    weights = np.array([nprim[e] for e in e_level[keep]])
    act = activities_df.loc[keep, time_cols].to_numpy(dtype=float) * weights[:, None]
    pair_act = (pd.DataFrame(act, index=activities_df.index[keep])
                  .groupby(level=["volume", "isotope"]).sum())
    if activity_threshold > 0:
        pair_act = pair_act.where(pair_act > activity_threshold, 0.0)
    pair_act = pair_act[(pair_act > 0).any(axis=1)]
    print(f"  {len(pair_act)} (volume, isotope) pairs with activity")

    # ---- Spectrum: sum over pairs of activity x spectrum per decay ------
    total = np.zeros((len(times_s), len(src.energy)))
    missing: list[str] = []
    rows: list[np.ndarray] = []
    specs: list[np.ndarray] = []

    def flush() -> None:
        if rows:
            np.add(total, np.array(rows).T @ np.array(specs), out=total)
            rows.clear()
            specs.clear()

    for (volume, isotope), a in zip(pair_act.index, pair_act.to_numpy()):
        spec = src.get(volume, isotope)
        if spec is None:
            missing.append(f"{volume}/{isotope}")
            continue
        rows.append(a)
        specs.append(spec)
        if len(rows) == 2000:
            flush()
    flush()

    print(f"  {len(pair_act) - len(missing)} pairs accumulated.")
    if missing:
        print(f"  [WARN] {len(missing)} spectrum(s) missing in {src.source}, "
              f"e.g. {', '.join(missing[:5])}")

    # ---- Package output -----------------------------------------------
    index = pd.Index(times_s, name="time_s")
    spectra_df = pd.DataFrame(total.astype(np.float32), index=index,
                              columns=[f"{e:.4f}" for e in src.energy])
    spectra_df.attrs.update({
        "timestep_s":         timestep,
        "duty":               duty,
        "beam_area_cm2":      beam_area,
        "nprim_per_step":     nprim,
        "in_belt_only":       in_belt_only,
        "activity_threshold": activity_threshold,
        "spectra_source":     src.source,
        "spectra_mode":       src.mode,
        "spectra_kind":       src.kind,
        "edges_keV":          src.edges.tolist(),
        "emin_keV":           emin,
        "emax_keV":           emax,
    })
    count_rate = pd.Series(band_integral(total, src.widths, mask), index=index,
                           name="count_rate_cps")
    count_rate.attrs.update({"spectra_source": src.source, "spectra_mode": src.mode,
                             "emin_keV": emin, "emax_keV": emax})
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

    header = "time_s   count_rate_cps"
    a = spectra_df.attrs
    if a.get("spectra_kind") == "npz" or a.get("emin_keV") is not None or a.get("emax_keV") is not None:
        header += (f"\nsource: {a.get('spectra_source')}\nmode: {a.get('spectra_mode')}"
                   f"\nband_keV: {a.get('emin_keV')} {a.get('emax_keV')}")
    np.savetxt(
        outdir / "count_rate.dat",
        np.column_stack([count_rate.index.to_numpy(), count_rate.to_numpy()]),
        header=header,
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
    band:       str             = "",
) -> None:
    """
    Log-log plot of total count rate [counts/s] vs decay time [s].
    Only time points with count rate > 0 are plotted.  band: energy band of the
    count rate, shown in the title.
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
    ax.set_title(title + (f" [{band}]" if band else ""))
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
    xscale:     str             = "log",
    band:       tuple[float | None, float | None] = (None, None),
) -> None:
    """
    Log-log plot of background spectra for every time point in spectra_df.

    Each spectrum is drawn as a step curve.  Curves are coloured by decay
    time using the "magma" colormap, with the colour proportional to
    log10(time_s) so that decades are evenly spaced in colour space.
    A colorbar shows the time axis.  Time points with all-zero spectra
    are silently skipped.  xscale: 'log' or 'linear' energy axis; band: (emin,
    emax) of the count rate, used as the x limits when given.
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

    ax.set_xscale(xscale)
    ax.set_yscale("log")
    lo, hi = band
    if lo is not None or hi is not None:
        ax.set_xlim(lo if lo is not None else (max(energies[0], 1.0) if xscale == "log"
                                               else energies[0]),
                    hi if hi is not None else energies[-1])
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
        description="Background spectra and count rate after one orbital time "
                    "step of irradiation at the orbit-mean proton spectrum.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("activities",   help="activities.pkl from compute_activities.py")
    p.add_argument("spenvis_file", help="SPENVIS AP9/AE9 output text file")
    p.add_argument("--spectra-dir",  default="result_spectra",  metavar="DIR",
                   help="Directory containing background spectrum .dat files "
                        "(fixed S-mode axis, 18-4112 keV).")
    p.add_argument("--spectra-file", default=None, metavar="NPZ",
                   help="spectra_<mode>.npz from build_spectra.py; replaces "
                        "--spectra-dir, and the energy axis comes from the file.")
    p.add_argument("--emin", type=float, default=None, metavar="KEV",
                   help="Lower edge of the energy band of the count rate [keV] "
                        "(default: whole axis).")
    p.add_argument("--emax", type=float, default=None, metavar="KEV",
                   help="Upper edge of the energy band of the count rate [keV].")
    p.add_argument("--xscale", choices=["log", "linear"], default="log",
                   help="Energy axis of the spectra plot.")
    p.add_argument("--energies", nargs="+", type=float, default=None, metavar="E",
                   help="Simulation energies to use [MeV] (default: all in activities).")
    p.add_argument("--threshold", type=float, default=0.0, metavar="BQ",
                   help="Skip (volume, isotope) activities after one step not above this [Bq]; "
                        "a speed-up only, it biases the long-term results.")
    p.add_argument("--R",         type=float, default=None, metavar="CM",
                   help="Source sphere radius (default: from the activities).")
    p.add_argument("--thetamax",  type=float, default=None, metavar="DEG",
                   help="Source cone half-angle (default: from the activities).")
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
        R=args.R,
        thetamax=args.thetamax,
        in_belt_only=args.in_belt_only,
        spectra_file=args.spectra_file,
        emin=args.emin,
        emax=args.emax,
    )

    print(f"\nTotal count rate summary:")
    nonzero_cr = count_rate[count_rate > 0]
    if len(nonzero_cr):
        print(f"  Peak    : {nonzero_cr.max():.4e} counts/s  "
              f"at t={nonzero_cr.idxmax():.2e} s")
        print(f"  Zero    : {count_rate[count_rate == 0].shape[0]} time points")
    else:
        print("  No contributions.")

    save_outputs(spectra_df, count_rate, outdir=args.outdir, fmt=args.fmt)
    plot_count_rate(count_rate, save_path=args.save_plot,
                    band=band_label(args.emin, args.emax))
    plot_spectra(spectra_df, save_path=args.save_spectra, xscale=args.xscale,
                 band=(args.emin, args.emax))


if __name__ == "__main__":
    main()
