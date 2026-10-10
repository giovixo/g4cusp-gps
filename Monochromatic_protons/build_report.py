#!/usr/bin/env python3
"""
Report on an activation-pipeline run: HTML page and PDF.

For each detector event class (default: scat_single, abs_single, compton) it
runs steps 6 and 8 on the step-4 spectra (or reuses their outputs; the step-7
histories are computed internally with activation_history), then collects

  - the instantaneous rate during the first day in orbit (60 s steps, belt passages marked);
  - the out-of-belt mean rate over a trailing 30-day window for the first 5 years;
  - the out-of-belt mean spectrum after one year;
  - the out-of-belt mean activity of every (volume, isotope) pair after 1 day, 1 month,
    1 year and 5 years in orbit, and its contribution to each event class;

and writes <outdir>/report_data.json, <outdir>/<name>.html (from report_template.html)
and <outdir>/<name>.pdf.

Usage (from Monochromatic_protons/):

    python build_report.py output/activities.pkl AP8MIN.AP8.output_mean_flux_550km_SSO.txt \\
        --spectra-dir result_postact --outdir report_run

Options: --modes, --reuse (skip steps 6 and 8 whose outputs exist), --name.
Steps 6 and 8 take a few minutes; the curves and statistics about two more.
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

from activation_history import compute_long_term_history     # noqa: E402
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
            if reuse and target.exists():
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


def summary(activities: Path, spenvis: Path, spectra_dir: Path, modes: list[str], ps: pd.DataFrame) -> dict:
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
        "date": dt.date.today().isoformat(), "commit": commit, "spenvis": spenvis.name,
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
    template = (HERE / "report_template.html").read_text()
    path.write_text(template.replace("__DATA__", json.dumps(data, separators=(",", ":"))))


def write_pdf(data: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    modes = list(data["summary"]["modes"])
    names = data["summary"]["modes"]
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    styles = ["-", "-", "--"]
    S, T = data["summary"], data["tables"]
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
                import textwrap
                for line in textwrap.wrap(content, 105):
                    fig.text(0.08, y, line, fontsize=9, va="top")
                    y -= 0.016
                y -= 0.01
            elif kind == "h":
                y -= 0.006
                fig.text(0.08, y, content, fontsize=11, weight="bold", va="top")
                y -= 0.024
            elif kind == "table":
                import textwrap
                head, rows = content
                # column widths from the longest entry, capped; long cells are wrapped
                char_w, max_chars = 0.0098, 85          # page fraction per character at 7.5 pt
                lens = [max(len(str(r[c])) for r in rows + [head]) + 2 for c in range(len(head))]
                widths = [max(n, 7) for n in lens]
                if sum(widths) > max_chars:             # shrink the widest columns and wrap
                    while sum(widths) > max_chars:
                        k = int(np.argmax(widths))
                        widths[k] -= 1
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

    with PdfPages(path) as pdf:
        sel = S["selection"]
        r1 = T["rates"]["1 year"]
        text_page(pdf, "Activation-induced background of the CUSP polarimeter", [
            ("p", f"Activation-pipeline run, {S['date']} (g4cusp-gps branch activation, commit {S['commit']}). "
                  f"Orbit: {S['spenvis']}. Background from the decay of nuclides produced by trapped protons in "
                  f"the CUSP mass model, computed with the three-step method of Campana et al. (2026). Rates are "
                  f"means over the out-of-belt part of the orbit unless stated otherwise."),
            ("h", "Out-of-belt mean rate after one year in orbit"),
            ("table", (["Event class", "Rate [counts/s]"], [[names[m], g(r1[m])] for m in modes]
                       + [["Payload activity", f"{g(r1['activity_Bq'])} Bq"]])),
            ("h", "Run summary"),
            ("table", (["Step", "Content"], [
                ["0-1", f"Nuclide production: {S['n_energies']} proton energies, {S['energies_MeV'][0]:g}-"
                        f"{S['energies_MeV'][1]:g} MeV, {S['nprim']:.0e} protons each"],
                ["2-3", "Decay chains and activities per volume and nuclide (Geant4 RadioactiveDecay6.1.2)"],
                ["4", f"Decays at rest: {S['n_decays']:.2e} decays, {S['n_pairs']} pairs "
                      f"({S['n_volumes']} volumes, {S['n_isotopes']} nuclides), {S['decays_min']:.0e}-"
                      f"{S['decays_max']:.0e} per pair"],
                ["4b", f"Selection: thresholds {sel['thr_scat']:g} keV plastic / {sel['thr_abs']:g} keV GAGG, "
                       f"window {sel['window']:g} ns, 0-{sel['emax']:g} keV ({sel['selection']})"],
                ["5-8", "Orbit history (60 s steps), 5-year history, out-of-belt averages"]])),
            ("h", "Out-of-belt mean rate by time in orbit"),
            ("table", (["Time in orbit", "Activity [Bq]"] + [f"{names[m]} [counts/s]" for m in modes],
                       [[t, g(v["activity_Bq"])] + [g(v[m]) for m in modes] for t, v in T["rates"].items()])),
        ])

        # first day
        fig, ax = plt.subplots(figsize=(A4[0], 4.6))
        for s, e in data["day"]["belt"]:
            ax.axvspan(s, e, color="#7a8496", alpha=0.18, lw=0)
        th = np.array(data["day"]["t_h"])
        for i, m in enumerate(modes):
            y = np.array(data["day"][m], dtype=float)
            ax.plot(th, np.where(y > 0, y, np.nan), styles[i], color=colors[i], lw=1.1, label=names[m])
        ax.set_yscale("log"); ax.set_ylim(bottom=1e-2); ax.set_xlim(0, 24)
        ax.set_xlabel("Time since start of mission [h]"); ax.set_ylabel("Rate [counts/s]")
        ax.set_title(pad=26, label="Activation-induced rate during the first day (shaded: in radiation belt)")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # 5 years
        fig, ax = plt.subplots(figsize=(A4[0], 4.6))
        td = np.array(data["long"]["t_d"])
        for i, m in enumerate(modes):
            ax.plot(td, data["long"][m], styles[i], color=colors[i], lw=1.6, label=names[m])
        ax.set_yscale("log"); ax.set_xlim(0, td[-1])
        ax.set_xlabel("Time since start of mission [d]"); ax.set_ylabel("Rate [counts/s]")
        ax.set_title(pad=26, label=f"Out-of-belt mean rate, trailing {LONG_WINDOW_D}-day window, first 5 years")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # spectra
        fig, ax = plt.subplots(figsize=(A4[0], 4.6))
        e = np.array(data["spec"]["e_keV"])
        for i, m in enumerate(modes):
            y = np.array(data["spec"][m], dtype=float)
            ax.step(e, np.where(y > 0, y, np.nan), where="mid", ls=styles[i], color=colors[i], lw=0.9, label=names[m])
        ax.set_yscale("log"); ax.set_xlim(0, e[-1] + data["spec"]["width_keV"] / 2); ax.set_ylim(bottom=1e-6)
        ax.set_xlabel("Deposited energy [keV]"); ax.set_ylabel("Rate [counts/s/keV]")
        ax.set_title(pad=26, label=f"Out-of-belt mean spectrum after one year ({data['spec']['width_keV']:g} keV bins)")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # tables
        for t in TIMES:
            text_page(pdf, f"Where the activity is: after {t} in orbit", [
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
            text_page(pdf, f"{names[m]}: background after one year", [
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
            ("p", "Orbit: the SPENVIS file repeats its period; the AP8MIN file stops at 400 MeV, so simulation "
                  "energies above it receive no protons."),
            ("p", "Coverage: the (volume, nuclide) pairs below the step-3 activity threshold (about 0.9% of the "
                  "out-of-belt activity in this run) were not simulated in step 4; rates are low by about that fraction."),
            ("p", "Geant4 data: decays follow Geant4's own decay data, including the known RadioactiveDecay6.1.2 "
                  "errors (Geant4 bug 2780 and related reports)."),
            ("p", "Statistics: the statistical error of the total rates is below 1e-3 (relative); single weak pairs "
                  "are noisier."),
        ])


# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("activities", type=Path, help="activities.pkl from compute_activities.py")
    p.add_argument("spenvis_file", type=Path, help="SPENVIS file used for the orbit")
    p.add_argument("--spectra-dir", type=Path, default=Path("result_postact"),
                   help="directory with spectra_<mode>.npz and the step-4 runs files")
    p.add_argument("--outdir", type=Path, default=Path("report_run"))
    p.add_argument("--modes", nargs="+", default=["scat_single", "abs_single", "compton"])
    p.add_argument("--reuse", action="store_true", help="reuse existing step 6 and 8 outputs in --outdir")
    p.add_argument("--name", default="cusp_activation_report", help="base name of the HTML and PDF files")
    args = p.parse_args(argv)

    activities, spenvis = args.activities.resolve(), args.spenvis_file.resolve()
    spectra_dir, outdir = args.spectra_dir.resolve(), args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("Steps 6 and 8 ...")
    run_steps(activities, spenvis, spectra_dir, outdir, args.modes, args.reuse)
    print("Curves ...")
    data = curves(outdir, spenvis, args.modes)
    print("Per-pair activities and rates ...")
    ps = pair_stats(activities, spenvis, spectra_dir, args.modes)
    ps.to_csv(outdir / "pair_stats.csv", index=False)
    data["tables"] = tables(ps, args.modes)
    data["summary"] = summary(activities, spenvis, spectra_dir, args.modes, ps)
    (outdir / "report_data.json").write_text(json.dumps(data, separators=(",", ":")))

    write_html(data, outdir / f"{args.name}.html")
    write_pdf(data, outdir / f"{args.name}.pdf")
    print(f"Report: {outdir / (args.name + '.html')}, {outdir / (args.name + '.pdf')}")


if __name__ == "__main__":
    main()
