#!/usr/bin/env python3
"""
Report on an activation-pipeline run: HTML page and PDF, for one or more orbits.

For each orbit (SPENVIS file) and each detector event class (default:
scat_single, abs_single, compton) it runs steps 6 and 8 on the step-4 spectra
(or reuses their outputs; the step-7 histories are computed internally with
activation_history), then collects

  - the instantaneous rate during the first day in orbit (60 s steps, belt passages marked);
  - the out-of-belt mean rate over a trailing 30-day window for the first 5 years;
  - the out-of-belt mean spectrum after one year;
  - the out-of-belt mean activity of every (volume, isotope) pair after 1 day, 1 month,
    1 year and 5 years in orbit, and its contribution to each event class;

and, to compare the orbits, their orbit-averaged proton spectra and the primaries
per simulation energy.  The step-0 to step-4 outputs do not depend on the orbit:
a new orbit or flux model needs only this script.  Writes <outdir>/report_data.json,
<outdir>/<name>.html (from report_template.html) and <outdir>/<name>.pdf; the step
6 and 8 outputs go to <outdir>/<label>/<mode>/.

Usage (from activation-pipeline/):

    python build_report.py output/activities.pkl AP8MIN.AP8.output_mean_flux_550km_SSO.txt \\
        AP9MEAN.AP9.output_mean_flux_550km_SSO.txt --spectra-dir result_postact --outdir report_run

Options: --labels (default: the file names up to the first dot, e.g. AP8MIN), --modes,
--reuse (skip steps 6 and 8 whose outputs exist), --name.  The first orbit is the
reference of the ratios.  About 5 minutes per orbit.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from accumulate_spectra import primaries_per_step, source_beam_area   # noqa: E402
from activation_history import (compute_long_term_history, count_rate_info,   # noqa: E402
                                flux_profile, read_flux_table)
from nuclides import parse_name                              # noqa: E402
from pair_spectra import load_pair_spectra, spectra_path     # noqa: E402
from run_postactivation import compute_pair_weights          # noqa: E402

DAY = 86400.0
YEAR = 365.25 * DAY
TIMES = {"1 day": DAY, "1 month": 30 * DAY, "1 year": YEAR, "5 years": 5 * YEAR}
NAMES = {"scat_single": "Plastic singles", "abs_single": "GAGG singles",
         "compton": "Compton coincidences", "any": "Any deposit"}
SPEC_REBIN = 4            # channels per displayed bin in the 1-year spectra
LONG_WINDOW_D = 30        # trailing window of the long-term curve [days]


# ---------------------------------------------------------------------------
# Steps 6-8
# ---------------------------------------------------------------------------

def run_steps(activities: Path, spenvis: Path, spectra_dir: Path, outdir: Path,
              modes: list[str], reuse: bool) -> None:
    """Steps 6 (count rate, spectra), 7 (1 day and 5 years) and 8 (1 year) per mode."""
    py = sys.executable
    for m in modes:
        d = outdir / m
        d.mkdir(parents=True, exist_ok=True)
        npz = spectra_path(spectra_dir, m)
        jobs = [
            (d / "count_rate.dat",
             [py, "accumulate_spectra.py", str(activities), str(spenvis), "--spectra-file", str(npz),
              "--outdir", str(d)], "step6"),
            (d / "avg_1y.dat",
             [py, "average_spectrum.py", str(d / "spectra.pkl"), str(spenvis), "--duration", "1y",
              "--activities", str(activities), "--spectra-file", str(npz), "--xscale", "linear",
              "--save-dat", str(d / "avg_1y.dat"), "--save-plot", str(d / "avg_1y.png")], "step8"),
        ]
        for target, cmd, tag in jobs:
            # outputs written before the flux time profile was recorded are redone
            if reuse and target.exists() and \
                    count_rate_info(d / "count_rate.dat").get("profile_emin") is not None:
                print(f"[{m}] {tag}: reusing {target}")
                continue
            print(f"[{m}] {tag}: {' '.join(cmd[1:3])} ...")
            with open(d / f"{tag}.log", "w") as log:
                subprocess.run(cmd, cwd=HERE, stdout=log, stderr=subprocess.STDOUT, check=True)


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def _sig(values, digits: int = 4) -> list[float]:
    return [float(f"{x:.{digits}g}") for x in values]


def curves(outdir: Path, spenvis: Path, modes: list[str]) -> dict:
    """First-day instantaneous rate, 5-year trailing out-of-belt mean, 1-year spectra."""
    out: dict = {}
    # first day
    day, belt, th = {}, None, None
    for m in modes:
        t, C, _, outb = compute_long_term_history(outdir / m / "count_rate.dat", spenvis, DAY, 3600.0)
        n = int(round(DAY / (t[1] - t[0]))) + 1
        day[m], belt, th = C[:n], ~outb[:n], t[:n] / 3600.0
    segs, start = [], None
    for i, b in enumerate(belt):
        if b and start is None:
            start = th[i]
        elif not b and start is not None:
            segs.append([round(start, 4), round(th[i], 4)])
            start = None
    if start is not None:
        segs.append([round(start, 4), round(th[-1], 4)])
    out["day"] = {"t_h": np.round(th, 4).tolist(), "belt": segs, **{m: _sig(day[m]) for m in modes}}

    # 5 years
    long = {}
    for m in modes:
        t, C, _, outb = compute_long_term_history(outdir / m / "count_rate.dat", spenvis, 5 * YEAR, DAY)
        per_day = int(round(DAY / (t[1] - t[0])))
        nd = int(5 * YEAR // DAY)
        Cs = C[:nd * per_day].reshape(nd, per_day)
        Os = outb[:nd * per_day].reshape(nd, per_day)
        num = np.convolve((Cs * Os).sum(1), np.ones(LONG_WINDOW_D))[:nd]
        den = np.convolve(Os.sum(1).astype(float), np.ones(LONG_WINDOW_D))[:nd]
        long[m] = num / den
    out["long"] = {"t_d": list(range(1, nd + 1)), **{m: _sig(long[m]) for m in modes}}

    # 1-year spectra
    spec, ec = {}, None
    for m in modes:
        a = np.loadtxt(outdir / m / "avg_1y.dat")
        n = len(a) // SPEC_REBIN
        spec[m] = a[:n * SPEC_REBIN, 1].reshape(n, SPEC_REBIN).mean(1)
        ec = a[:n * SPEC_REBIN, 0].reshape(n, SPEC_REBIN).mean(1)
        widths = np.diff(a[:2, 0])[0] * SPEC_REBIN
    out["spec"] = {"e_keV": np.round(ec, 2).tolist(), "width_keV": float(widths), **{m: _sig(spec[m]) for m in modes}}
    return out


def flux_data(activities: Path, spenvis: Path) -> dict:
    """Orbit-averaged integral spectrum, primaries per simulation energy, time-profile statistics."""
    levels, flux, _ = read_flux_table(spenvis)
    acts_attrs = pd.read_pickle(activities).attrs
    energies = sorted(float(e) for e in acts_attrs.get("nprim", {}))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # bands above the last level: reported below
        nprim, _, _ = primaries_per_step(spenvis, energies, source_beam_area(acts_attrs))
    total, prof = flux[:, 0], flux_profile(spenvis, energies[0])
    return {
        "levels_MeV": levels.tolist(), "F_mean": _sig(flux.mean(axis=0)),
        "nprim": {f"{e:g}": float(f"{nprim[e]:.4g}") for e in energies},
        "emax_MeV": float(levels[-1]), "profile_emin_MeV": energies[0],
        "peak_to_mean_total": round(float(total.max() / total.mean()), 1),
        "peak_to_mean_profile": round(float(prof.max() / prof.mean()), 1),
        "in_belt_fraction": round(float((total > 0).mean()), 4),
    }


def pair_stats(activities: Path, spenvis: Path, spectra_dir: Path, modes: list[str]) -> pd.DataFrame:
    """Out-of-belt mean activity of each pair at the TIMES, and its rate in each mode."""
    libs = {m: load_pair_spectra(spectra_path(spectra_dir, m)) for m in modes}
    rows = []
    for label, T in TIMES.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = compute_pair_weights(activities, spenvis, duration_s=T)
        w["time"] = label
        idx = pd.MultiIndex.from_arrays([w.Volume, w.Isotope])
        for m, lib in libs.items():
            eff = np.asarray(lib.counts.sum(axis=1)).ravel() / lib.ndecays
            e = pd.Series(eff, index=pd.MultiIndex.from_arrays([lib.volumes, lib.isotopes]))
            w[f"rate_{m}"] = w["Weight_Bq"].to_numpy() * e.reindex(idx).fillna(0.0).to_numpy()
        rows.append(w)
        print(f"  {label}: activity {w.Weight_Bq.sum():.4g} Bq, "
              + ", ".join(f"{m} {w[f'rate_{m}'].sum():.4g} counts/s" for m in modes))
    return pd.concat(rows, ignore_index=True)


def volume_group(volume: str) -> str:
    rules = [(r"^PV-Scatterer_\d{3}$", "Plastic scatterers"), (r"^PV-Absorber_\d{3}$", "GAGG absorbers"),
             (r"Wrapping", "Scintillator wrappings"), (r"PCB|MAPMT", "PCBs and photosensors"),
             (r"[Ff]ilter|_Collimator", "Filters and collimators"), (r"CUSP_GEANT4", "World volume")]
    for pattern, name in rules:
        if re.search(pattern, volume):
            return name
    return "Structure (frames, lids, panels, truss)"


def _half_lives() -> dict:
    data = Path(os.environ.get("GEANT4_DATA_DIR", "")) / "G4ENSDFSTATE3.0" / "ENSDFSTATE.dat"
    table = {}
    if data.exists():
        for line in open(data):
            t = line.split()
            if len(t) >= 5 and float(t[4]) > 0:
                table[(int(t[0]), int(t[1]), round(float(t[2]), 1), t[3].lstrip("+-").upper())] = \
                    float(t[4]) * 1e-9 * math.log(2)
    return table


def fmt_time(s: float | None) -> str:
    if s is None:
        return "–"
    for unit, sec in (("y", YEAR), ("d", DAY), ("h", 3600.0), ("min", 60.0)):
        if s >= sec:
            return f"{s / sec:.3g} {unit}"
    return f"{s:.3g} s"


def tables(ps: pd.DataFrame, modes: list[str]) -> dict:
    hl_table = _half_lives()

    def half_life(iso: str) -> str:
        try:
            Z, A, E, flb = parse_name(iso)[:4]
        except Exception:
            return "–"
        for key in [(Z, A, round(E, 1), (flb or "").upper()), (Z, A, round(E, 1), "")]:
            if key in hl_table:
                return fmt_time(hl_table[key])
        return "–"

    ps = ps.assign(Group=ps.Volume.map(volume_group))
    out: dict = {"rates": {}, "volumes": {}, "isotopes": {}, "groups": {}, "contrib": {}}
    for T in TIMES:
        d = ps[ps.time == T]
        tot = d.Weight_Bq.sum()
        out["rates"][T] = {"activity_Bq": float(f"{tot:.4g}"), **{m: float(f"{d[f'rate_{m}'].sum():.4g}") for m in modes}}
        v = d.groupby("Volume").Weight_Bq.sum().nlargest(10)
        out["volumes"][T] = [[k, float(f"{x:.4g}"), round(100 * x / tot, 1)] for k, x in v.items()]
        i = d.groupby("Isotope").Weight_Bq.sum().nlargest(10)
        out["isotopes"][T] = [[k, half_life(k), float(f"{x:.4g}"), round(100 * x / tot, 1)] for k, x in i.items()]
        g = d.groupby("Group")[["Weight_Bq"] + [f"rate_{m}" for m in modes]].sum()
        out["groups"][T] = {name: {"act": round(100 * r.Weight_Bq / tot, 1),
                                   **{m: round(100 * r[f"rate_{m}"] / max(d[f'rate_{m}'].sum(), 1e-300), 1)
                                      for m in modes}}
                            for name, r in g.iterrows() if name != "World volume"}
    d = ps[ps.time == "1 year"]
    for m in modes:
        tot = d[f"rate_{m}"].sum()
        p = d.nlargest(10, f"rate_{m}")
        out["contrib"][m] = [[r.Volume, r.Isotope, float(f"{r[f'rate_{m}']:.3g}"),
                              round(100 * r[f"rate_{m}"] / tot, 1)] for _, r in p.iterrows()]
        iso = d.groupby("Isotope")[f"rate_{m}"].sum().nlargest(8)
        out["contrib"][m + "_iso"] = [[k, half_life(k), float(f"{x:.3g}"), round(100 * x / tot, 1)]
                                      for k, x in iso.items()]
    return out


def summary(activities: Path, spectra_dir: Path, modes: list[str], orbits: list[dict]) -> dict:
    attrs = pd.read_pickle(activities).attrs
    nprim = attrs.get("nprim", {})
    energies = sorted(float(e) for e in nprim)
    meta = load_pair_spectra(spectra_path(spectra_dir, modes[0])).meta
    runs = pd.concat([pd.read_csv(f) for f in sorted(spectra_dir.glob("batch_*_runs.csv"))])
    try:
        commit = subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                                         text=True).strip()
    except Exception:
        commit = "unknown"
    return {
        "date": dt.date.today().isoformat(), "commit": commit, "orbits": orbits,
        "energies_MeV": [energies[0], energies[-1]] if energies else None, "n_energies": len(energies),
        "nprim": int(np.median(list(nprim.values()))) if nprim else None,
        "source": attrs.get("source", {}),
        "n_pairs": int(len(runs)), "n_decays": int(runs.NDecays.sum()),
        "n_volumes": int(runs.Volume.nunique()), "n_isotopes": int(runs.Isotope.nunique()),
        "decays_min": int(runs.NDecays.min()), "decays_max": int(runs.NDecays.max()),
        "selection": {k: meta.get(k) for k in ("thr_scat", "thr_abs", "window", "fwhm_scat", "fwhm_abs",
                                               "emax", "binwidth", "selection")},
        "modes": {m: NAMES.get(m, m) for m in modes},
    }


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def write_html(data: dict, path: Path) -> None:
    template = (HERE / "report_template.html").read_text(encoding="utf-8")
    path.write_text(template.replace("__DATA__", json.dumps(data, separators=(",", ":"))), encoding="utf-8")


def write_pdf(data: dict, path: Path) -> None:
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    S = data["summary"]
    modes, names = list(S["modes"]), S["modes"]
    orbits = data["orbits"]
    labels = [o["label"] for o in orbits]
    ref = labels[0]
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]          # event classes
    ocolors = ["#3b4250", "#8a4fd8", "#c2185b"]         # orbits (flux plot)
    ostyles = ["-", "--", ":"]                           # orbits (rate plots)
    plt.style.use("default")          # independent of the user's matplotlibrc
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.labelsize": 10,
                         "axes.titlesize": 11, "axes.titleweight": "bold", "legend.fontsize": 9,
                         "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "axes.grid": True,
                         "grid.color": "#e3e6eb", "grid.linewidth": 0.6, "axes.edgecolor": "#9aa1ab",
                         "text.usetex": False, "mathtext.default": "regular"})
    A4 = (8.27, 11.69)

    def text_page(pdf, title, blocks):
        fig = plt.figure(figsize=A4)
        fig.text(0.08, 0.95, title, fontsize=15, weight="bold", va="top")
        y = 0.91
        for kind, content in blocks:
            if kind == "p":
                for line in textwrap.wrap(content, 105):
                    fig.text(0.08, y, line, fontsize=9, va="top")
                    y -= 0.016
                y -= 0.01
            elif kind == "h":
                y -= 0.006
                fig.text(0.08, y, content, fontsize=11, weight="bold", va="top")
                y -= 0.024
            elif kind == "table":
                head, rows = content
                # column widths from the longest entry, capped; long cells are wrapped
                char_w, max_chars = 0.0098, 85          # page fraction per character at 7.5 pt
                lens = [max(len(str(r[c])) for r in rows + [head]) + 2 for c in range(len(head))]
                widths = [max(n, 7) for n in lens]
                while sum(widths) > max_chars:          # shrink the widest columns and wrap
                    widths[int(np.argmax(widths))] -= 1
                total = sum(widths)
                frac = [w / total for w in widths]
                width = char_w * total
                wrap = lambda x, c: textwrap.fill(str(x), max(widths[c] - 2, 4))
                head = [wrap(x, c) for c, x in enumerate(head)]
                wrapped = [[wrap(x, c) for c, x in enumerate(r)] for r in rows]
                nlines = [max(cell.count("\n") + 1 for cell in r) for r in wrapped]
                nhead = max(x.count("\n") + 1 for x in head)
                line_h = 0.0145
                h = line_h * (nhead + 0.25 + sum(n + 0.25 for n in nlines))
                ax = fig.add_axes([0.08, y - h, width, h])
                ax.axis("off")
                tb = ax.table(cellText=wrapped, colLabels=head, colWidths=frac, loc="upper left",
                              cellLoc="left", colLoc="left", bbox=[0, 0, 1, 1])
                tb.auto_set_font_size(False)
                tb.set_fontsize(7.5)
                heights = [nhead + 0.25] + [n + 0.25 for n in nlines]
                for (r, c), cell in tb.get_celld().items():
                    cell.set_edgecolor("#dfe3e9")
                    cell.set_linewidth(0.5)
                    cell.set_height(heights[r] / sum(heights))
                    if r == 0:
                        cell.set_text_props(weight="bold", color="#4a515c")
                y -= h + 0.025
        pdf.savefig(fig)
        plt.close(fig)

    def g(x):
        return f"{x:.3g}"

    def rate_rows(t):
        """One row per orbit at time t, and the ratios to the reference orbit."""
        rows = []
        for o in orbits:
            v = o["tables"]["rates"][t]
            rows.append([o["label"], g(v["activity_Bq"])] + [g(v[m]) for m in modes])
        r0 = orbits[0]["tables"]["rates"][t]
        for o in orbits[1:]:
            v = o["tables"]["rates"][t]
            rows.append([f"{o['label']} / {ref}", f"{v['activity_Bq'] / r0['activity_Bq']:.3f}"]
                        + [f"{v[m] / r0[m]:.3f}" for m in modes])
        return rows

    sel = S["selection"]
    flux = data["flux"]
    e0 = S["energies_MeV"]
    with PdfPages(path) as pdf:
        text_page(pdf, "Activation-induced background of the CUSP polarimeter", [
            ("p", f"Activation-pipeline run, {S['date']} (g4cusp-gps branch activation, commit {S['commit']}). "
                  f"Orbits: " + "; ".join(f"{o['label']} ({o['spenvis']})" for o in S["orbits"]) + ". "
                  f"Background from the decay of nuclides produced by trapped protons in the CUSP mass model, "
                  f"computed with the three-step method of Campana et al. (2026). Rates are means over the "
                  f"out-of-belt part of the orbit unless stated otherwise. The simulation (steps 0-4) is the same "
                  f"for all orbits: only the proton flux changes."),
            ("h", "Out-of-belt mean rate after one year in orbit"),
            ("table", (["Orbit", "Activity [Bq]"] + [f"{names[m]} [counts/s]" for m in modes], rate_rows("1 year"))),
            ("h", "Run summary"),
            ("table", (["Step", "Content"], [
                ["0-1", f"Nuclide production: {S['n_energies']} proton energies, {e0[0]:g}-{e0[1]:g} MeV, "
                        f"{S['nprim']:.0e} protons each"],
                ["2-3", "Decay chains and activities per volume and nuclide (Geant4 RadioactiveDecay6.1.2)"],
                ["4", f"Decays at rest: {S['n_decays']:.2e} decays, {S['n_pairs']} pairs "
                      f"({S['n_volumes']} volumes, {S['n_isotopes']} nuclides), {S['decays_min']:.0e}-"
                      f"{S['decays_max']:.0e} per pair"],
                ["4b", f"Selection: thresholds {sel['thr_scat']:g} keV plastic / {sel['thr_abs']:g} keV GAGG, "
                       f"window {sel['window']:g} ns, 0-{sel['emax']:g} keV ({sel['selection']})"],
                ["5-8", "Per orbit: history at 60 s steps, 5-year history, out-of-belt averages"]])),
            ("h", "Out-of-belt mean rate by time in orbit"),
            ("table", (["Time in orbit", "Orbit", "Activity [Bq]"] + [f"{names[m]} [counts/s]" for m in modes],
                       [[t] + r for t in TIMES for r in rate_rows(t)])),
        ])

        # proton spectra
        fig, (ax, axr) = plt.subplots(2, 1, figsize=(A4[0], 7.4), sharex=True,
                                      gridspec_kw={"height_ratios": [2.2, 1]})
        diff = {}
        for k, lab in enumerate(labels):
            L, F = np.array(flux[lab]["levels_MeV"]), np.array(flux[lab]["F_mean"], float)
            ec, dfl = np.sqrt(L[:-1] * L[1:]), (F[:-1] - F[1:]) / np.diff(L)
            ok = dfl > 0
            diff[lab] = dict(zip(np.round(ec, 6), dfl))
            ax.plot(ec[ok], dfl[ok], "o-", color=ocolors[k], ms=3.5, lw=1.4, label=lab)
        ax.axvspan(e0[0], e0[1], color="#7a8496", alpha=0.10, lw=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylabel("Differential flux [p/cm²/s/MeV]")
        ax.set_title(pad=26, label="Orbit-averaged trapped-proton spectrum (shaded: simulated energies)")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, frameon=False)
        for k, lab in enumerate(labels[1:], 1):
            common = sorted(set(diff[ref]) & set(diff[lab]))
            x = np.array(common)
            y = np.array([diff[lab][c] / diff[ref][c] for c in common])
            axr.plot(x, y, "o-", color=ocolors[k], ms=3.5, lw=1.4, label=f"{lab} / {ref}")
        axr.axhline(1, color="#9aa1ab", lw=0.8)
        axr.axvspan(e0[0], e0[1], color="#7a8496", alpha=0.10, lw=0)
        axr.set_xscale("log")
        axr.set_xlabel("Proton energy [MeV]"); axr.set_ylabel(f"Ratio to {ref}")
        axr.legend(loc="upper left", frameon=False)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        Es = list(flux[ref]["nprim"])
        text_page(pdf, "Proton flux of the orbits", [
            ("p", "Orbit-averaged integral flux F(>E) of the SPENVIS files, and the primaries per 60 s step "
                  "assigned to each simulation energy (band between the geometric midpoints to its neighbours). "
                  "The time profile of the flux along the orbit is F(>E_min), with E_min the lowest simulation "
                  "energy; the belt passages are the steps with a non-zero total flux."),
            ("h", "Integral flux F(>E) [p/cm²/s]"),
            ("table", (["E [MeV]"] + labels + [f"{l} / {ref}" for l in labels[1:]],
                       [[f"{E:g}"] + [g(_fgt(flux[l], E)) if _fgt(flux[l], E) is not None else "–" for l in labels]
                        + [_ratio(_fgt(flux[l], E), _fgt(flux[ref], E)) for l in labels[1:]]
                        for E in (1, 6, 10, 20, 30, 60, 100, 200, 400, 700)])),
            ("h", "Primaries per 60 s step"),
            ("table", (["Simulation energy [MeV]"] + labels + [f"{l} / {ref}" for l in labels[1:]],
                       [[e] + [g(flux[l]["nprim"][e]) for l in labels]
                        + [_ratio(flux[l]["nprim"][e], flux[ref]["nprim"][e]) for l in labels[1:]] for e in Es]
                       + [["Total"] + [g(sum(flux[l]["nprim"].values())) for l in labels]
                          + [_ratio(sum(flux[l]["nprim"].values()), sum(flux[ref]["nprim"].values()))
                             for l in labels[1:]]])),
            ("h", "Time profile"),
            ("table", (["Orbit", "Last level [MeV]", "In-belt steps", "Peak/mean, total flux",
                        f"Peak/mean, F(>{flux[ref]['profile_emin_MeV']:g} MeV)"],
                       [[l, f"{flux[l]['emax_MeV']:g}", f"{100 * flux[l]['in_belt_fraction']:.1f}%",
                         g(flux[l]["peak_to_mean_total"]), g(flux[l]["peak_to_mean_profile"])] for l in labels])),
        ])

        # first day, one panel per orbit
        fig, axes = plt.subplots(len(orbits), 1, figsize=(A4[0], 3.6 * len(orbits) + 0.6), squeeze=False)
        for ax, o in zip(axes[:, 0], orbits):
            for s_, e_ in o["day"]["belt"]:
                ax.axvspan(s_, e_, color="#7a8496", alpha=0.18, lw=0)
            th = np.array(o["day"]["t_h"])
            for i, m in enumerate(modes):
                y = np.array(o["day"][m], dtype=float)
                ax.plot(th, np.where(y > 0, y, np.nan), color=colors[i], lw=1.1, label=names[m])
            ax.set_yscale("log"); ax.set_ylim(1e-2, None); ax.set_xlim(0, 24)
            ax.set_ylabel("Rate [counts/s]")
            ax.set_title(f"{o['label']}: rate during the first day (shaded: in radiation belt)", pad=22)
            ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, fontsize=8)
        axes[-1, 0].set_xlabel("Time since start of mission [h]")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # 5 years: colour = event class, line style = orbit
        fig, ax = plt.subplots(figsize=(A4[0], 4.8))
        for k, o in enumerate(orbits):
            td = np.array(o["long"]["t_d"])
            for i, m in enumerate(modes):
                ax.plot(td, o["long"][m], ostyles[k], color=colors[i], lw=1.5,
                        label=f"{names[m]}, {o['label']}")
        ax.set_yscale("log"); ax.set_xlim(0, td[-1])
        ax.set_xlabel("Time since start of mission [d]"); ax.set_ylabel("Rate [counts/s]")
        ax.set_title(pad=40, label=f"Out-of-belt mean rate, trailing {LONG_WINDOW_D}-day window, first 5 years")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(modes), frameon=False, fontsize=8)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # spectra: one panel per event class
        fig, axes = plt.subplots(len(modes), 1, figsize=(A4[0], 3.3 * len(modes) + 0.6), sharex=True, squeeze=False)
        for ax, (i, m) in zip(axes[:, 0], enumerate(modes)):
            for k, o in enumerate(orbits):
                e = np.array(o["spec"]["e_keV"]); y = np.array(o["spec"][m], dtype=float)
                ax.step(e, np.where(y > 0, y, np.nan), where="mid", ls=ostyles[k], color=colors[i], lw=0.9,
                        label=o["label"])
            ax.set_yscale("log"); ax.set_ylim(1e-6, None); ax.set_xlim(0, e[-1] + o["spec"]["width_keV"] / 2)
            ax.set_ylabel("Rate [counts/s/keV]")
            ax.set_title(f"{names[m]}: out-of-belt mean spectrum after one year", pad=6)
            ax.legend(loc="upper right", frameon=False)
        axes[-1, 0].set_xlabel("Deposited energy [keV]")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # tables, per orbit
        for o in orbits:
            T = o["tables"]
            for t in TIMES:
                text_page(pdf, f"{o['label']}: where the activity is after {t} in orbit", [
                    ("h", "Most active volumes"),
                    ("table", (["Volume", "Activity [Bq]", "Share"], [[v, g(a), f"{s}%"] for v, a, s in T["volumes"][t]])),
                    ("h", "Most active nuclides"),
                    ("table", (["Nuclide", "T1/2", "Activity [Bq]", "Share"], [[i, h, g(a), f"{s}%"] for i, h, a, s in T["isotopes"][t]])),
                    ("h", "Share by volume group (activity, and each event class)"),
                    ("table", (["Volume group", "Activity"] + [names[m] for m in modes],
                               [[k, f"{v['act']}%"] + [f"{v[m]}%" for m in modes]
                                for k, v in sorted(T["groups"][t].items(), key=lambda kv: -kv[1]["act"])])),
                ])
            for m in modes:
                text_page(pdf, f"{o['label']}, {names[m]}: background after one year", [
                    ("h", "Leading (volume, nuclide) pairs"),
                    ("table", (["Volume", "Nuclide", "Rate [counts/s]", "Share"],
                               [[v, i, g(r), f"{s}%"] for v, i, r, s in T["contrib"][m]])),
                    ("h", "Leading nuclides (all volumes)"),
                    ("table", (["Nuclide", "T1/2", "Rate [counts/s]", "Share"],
                               [[i, h, g(r), f"{s}%"] for i, h, r, s in T["contrib"][m + "_iso"]]))])
        text_page(pdf, "Method and limits", [
            ("p", "Event classes (provisional selection): a scintillator triggers above the plastic or GAGG threshold; "
                  "deposits later than the coincidence window after the first trigger are dropped. Plastic singles: "
                  "one plastic bar and no GAGG crystal. GAGG singles: one GAGG crystal and no plastic bar. Compton "
                  "coincidences: exactly one plastic bar and one GAGG crystal, energies summed; no geometric or "
                  "kinematic cuts yet."),
            ("p", "Detector response: true energy deposits; no energy resolution, light yield or quenching; no dead "
                  "time or pile-up. The prompt background during belt passages is not part of this calculation."),
            ("p", orbit_note(data)),
            ("p", "Coverage: the (volume, nuclide) pairs below the step-3 activity threshold (about 0.9% of the "
                  "out-of-belt activity) were not simulated in step 4; rates are low by about that fraction."),
            ("p", "Geant4 data: decays follow Geant4's own decay data, including the known RadioactiveDecay6.1.2 "
                  "errors (Geant4 bug 2780 and related reports)."),
            ("p", "Statistics: the statistical error of the total rates is below 1e-3 (relative); single weak pairs "
                  "are noisier. The same step-4 spectra serve all orbits, so the ratios between orbits carry "
                  "almost no statistical error."),
        ])


def _fgt(f: dict, E: float) -> float | None:
    """Orbit-mean F(>E) at a tabulated level, None if E is not one."""
    L = f["levels_MeV"]
    return f["F_mean"][L.index(E)] if E in L else None


def _ratio(a, b) -> str:
    return "–" if a is None or b in (None, 0) else f"{a / b:.3f}"


def orbit_note(data: dict) -> str:
    flux, e0 = data["flux"], data["summary"]["energies_MeV"]
    parts = []
    for lab, f in flux.items():
        if f["emax_MeV"] < e0[1]:
            parts.append(f"{lab} stops at {f['emax_MeV']:g} MeV, so protons above it are missing and the "
                         f"{e0[1]:g} MeV simulation energy receives none")
        else:
            parts.append(f"{lab} reaches {f['emax_MeV']:g} MeV; protons above {e0[1]:g} MeV (the highest "
                         f"simulation energy) are not included")
    return ("Orbit: each SPENVIS file repeats its period. The flux time profile is F(>"
            f"{flux[next(iter(flux))]['profile_emin_MeV']:g} MeV), the activating protons; the belt passages "
            "are the steps with a non-zero total flux. " + "; ".join(parts) + ".")


# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("activities", type=Path, help="activities.pkl from compute_activities.py")
    p.add_argument("spenvis_files", type=Path, nargs="+",
                   help="SPENVIS files, one per orbit or flux model (the first is the reference)")
    p.add_argument("--labels", nargs="+", default=None,
                   help="orbit labels (default: file names up to the first dot)")
    p.add_argument("--spectra-dir", type=Path, default=Path("result_postact"),
                   help="directory with spectra_<mode>.npz and the step-4 runs files")
    p.add_argument("--outdir", type=Path, default=Path("report_run"))
    p.add_argument("--modes", nargs="+", default=["scat_single", "abs_single", "compton"])
    p.add_argument("--reuse", action="store_true", help="reuse existing step 6 and 8 outputs in --outdir")
    p.add_argument("--name", default="cusp_activation_report", help="base name of the HTML and PDF files")
    args = p.parse_args(argv)

    activities = args.activities.resolve()
    spenvis_files = [f.resolve() for f in args.spenvis_files]
    labels = args.labels or [f.name.split(".")[0] for f in spenvis_files]
    if len(labels) != len(spenvis_files) or len(set(labels)) != len(labels):
        sys.exit("ERROR: give one distinct label per SPENVIS file")
    spectra_dir, outdir = args.spectra_dir.resolve(), args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    orbits, flux, all_ps = [], {}, []
    for label, spenvis in zip(labels, spenvis_files):
        odir = outdir / label
        print(f"=== {label}: {spenvis.name}")
        print("Steps 6 and 8 ...")
        run_steps(activities, spenvis, spectra_dir, odir, args.modes, args.reuse)
        print("Curves ...")
        o = {"label": label, "spenvis": spenvis.name, **curves(odir, spenvis, args.modes)}
        print("Per-pair activities and rates ...")
        ps = pair_stats(activities, spenvis, spectra_dir, args.modes)
        all_ps.append(ps.assign(orbit=label))
        o["tables"] = tables(ps, args.modes)
        orbits.append(o)
        flux[label] = flux_data(activities, spenvis)
    pd.concat(all_ps, ignore_index=True).to_csv(outdir / "pair_stats.csv", index=False)
    data = {"summary": summary(activities, spectra_dir, args.modes,
                               [{"label": o["label"], "spenvis": o["spenvis"]} for o in orbits]),
            "flux": flux, "orbits": orbits}
    data["summary"]["orbit_note"] = orbit_note(data)
    (outdir / "report_data.json").write_text(json.dumps(data, separators=(",", ":")))

    write_html(data, outdir / f"{args.name}.html")
    write_pdf(data, outdir / f"{args.name}.pdf")
    print(f"Report: {outdir / (args.name + '.html')}, {outdir / (args.name + '.pdf')}")


if __name__ == "__main__":
    main()
