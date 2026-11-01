#!/usr/bin/env python3
"""
Quick consistency check of a Geant4 RadioactiveDecay data library.

  1. Detail lines of non-IT decays (intensity > 0) with Q <= 0 (5th column, keV).
     G4VRadioactiveDecay::LoadDecayTable passes this column to the decay
     channel as the Q-value; with Q <= 0 no beta particle is emitted.
     The fraction of decays affected is computed the way LoadDecayTable does
     it: the detail intensities of each mode are rescaled to the mode's
     summary fraction.
  2. 'P'-line half-lives differing from G4ENSDFSTATE by more than 1 %
     (by default).  LoadDecayTable ignores them, as Geant4 takes the
     lifetimes from G4ENSDFSTATE, so they are informational.  Levels are
     paired by excitation energy (within 0.01 keV) and floating-level flag;
     P lines without a matching ENSDFSTATE level are not reported.
  3. Half-lives in the comment header of each file (first line, e.g.
     "# 90Y (  64.05 H  )") differing from the ENSDFSTATE ground state (E = 0,
     no floating flag) by more than the same tolerance.  Informational too;
     files whose header cannot be parsed, or whose ground state has no
     finite ENSDFSTATE half-life, are skipped.

Standard library only.
Usage:
    python3 check_rdm_quick.py [RadioactiveDecay dir] [ENSDFSTATE.dat] [--tol 0.01]
Defaults: $G4RADIOACTIVEDATA and $G4ENSDFSTATEDATA/ENSDFSTATE.dat (set by
geant4.sh / geant4.csh), else under $GEANT4_DATA_DIR.
Exit status: 0 if all three lists are empty, 1 otherwise.
"""
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

TOL_E = 1e-2   # keV, level-energy tolerance


def usage(msg):
    sys.exit(f"{msg}\n{__doc__}")


argv = sys.argv[1:]
tol_t = 0.01
if "--tol" in argv:
    i = argv.index("--tol")
    try:
        tol_t = float(argv[i + 1])
    except (IndexError, ValueError):
        usage("--tol needs a number")
    del argv[i:i + 2]
if any(a.startswith("-") for a in argv):
    usage(f"unknown option: {[a for a in argv if a.startswith('-')][0]}")

data = Path(os.environ.get("GEANT4_DATA_DIR", "."))


def default_dir(env, prefix):
    if os.environ.get(env):
        return Path(os.environ[env])
    found = sorted(data.glob(prefix + "*"))
    return found[-1] if found else None


rdm = Path(argv[0]) if len(argv) > 0 else default_dir("G4RADIOACTIVEDATA", "RadioactiveDecay")
if len(argv) > 1:
    ensdf_file = Path(argv[1])
else:
    d = default_dir("G4ENSDFSTATEDATA", "G4ENSDFSTATE")
    ensdf_file = d / "ENSDFSTATE.dat" if d else None
if rdm is None or not rdm.is_dir():
    usage(f"RadioactiveDecay directory not found: {rdm}")
if ensdf_file is None or not ensdf_file.is_file():
    usage(f"ENSDFSTATE.dat not found: {ensdf_file}")

SYM = ("n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge "
       "As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm "
       "Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U "
       "Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


def name(Z, A, E=0.0, flb=""):
    s = f"{SYM[Z]}{A}"
    return s + (f"[{E:g}{flb}]" if E > 0 or flb else "")


def fmt_t(t):
    for unit, sec in (("y", 3.15576e7), ("d", 86400), ("h", 3600), ("min", 60)):
        if t >= sec:
            return f"{t / sec:.4g} {unit}"
    return f"{t:.4g} s"


# ---- ENSDFSTATE: (Z, A) -> [(E keV, flb, T1/2 s or -1)] ----------------------
# columns: Z, A, E (keV), flag, mean life (ns; <= 0 for stable / unknown), ...
ensdf = defaultdict(list)
for line in open(ensdf_file):
    t = line.split()
    if len(t) >= 5:
        tau_ns = float(t[4])
        ensdf[(int(t[0]), int(t[1]))].append(
            (float(t[2]), t[3].lstrip("+-").upper(), tau_ns * 1e-9 * math.log(2) if tau_ns > 0 else -1))

# ---- RadioactiveDecay: (Z, A) -> [level dict] in file order --------------------
# P lines:       P, E (keV), flag, T1/2 (s)
# summary lines: mode, 0, fraction (3 columns)
# detail lines:  mode, daughter E (keV), flag, intensity (%), Q (keV) [, beta type]
lib = {}
for f in rdm.glob("z*.a*"):
    m = re.fullmatch(r"z(\d+)\.a(\d+)", f.name)
    if not m:
        continue
    levels = []
    for line in f.read_text().splitlines():
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "P":
            levels.append({"E": float(t[1]), "flb": t[2].lstrip("+-").upper(), "T": float(t[3]),
                           "summary": {}, "detail": []})
        elif levels and len(t) == 3:
            levels[-1]["summary"][t[0]] = float(t[2])
        elif levels and len(t) >= 5:
            levels[-1]["detail"].append((t[0], float(t[1]), float(t[3]), float(t[4])))
    lib[(int(m.group(1)), int(m.group(2)))] = levels

print(f"RadioactiveDecay: {rdm}\nENSDFSTATE:       {ensdf_file}\n")

print("1. Levels with non-IT decay branches having Q <= 0 (intensity > 0)")
n_q = 0
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        bad = [d for d in lv["detail"] if d[0] != "IT" and d[2] > 0 and d[3] <= 0]
        if bad:
            n_q += 1
            msum = defaultdict(float)
            for d in lv["detail"]:
                msum[d[0]] += d[2]
            frac = 100 * sum(lv["summary"].get(d[0], 0) * d[2] / msum[d[0]] for d in bad)
            modes = sorted({d[0] for d in bad})
            qs = ", ".join(f"{d[3]:g}" for d in bad[:3])
            print(f"   {name(Z, A, lv['E'], lv['flb']):16s} {len(bad):3d} branches, {frac:8.3f} % "
                  f"of decays, {'/'.join(modes)}; Q = {qs}{', ...' if len(bad) > 3 else ''}")
print(f"   total: {n_q} levels")

print(f"\n2. 'P'-line half-lives differing from ENSDFSTATE by more than {100 * tol_t:g} %")
n_t = 0
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        if lv["T"] <= 0:
            continue
        for E, flb, hl in ensdf.get((Z, A), []):
            if abs(E - lv["E"]) < TOL_E and flb == lv["flb"]:
                if hl > 0 and abs(lv["T"] / hl - 1) > tol_t:
                    n_t += 1
                    print(f"   {name(Z, A, lv['E'], lv['flb']):16s} P line {fmt_t(lv['T']):11s} "
                          f"ENSDFSTATE {fmt_t(hl):11s} ratio {lv['T'] / hl:.3g}")
                break
print(f"   total: {n_t} levels")

print(f"\n3. Comment-header half-lives differing from the ENSDFSTATE ground state by more than "
      f"{100 * tol_t:g} %")
UNIT = {"PS": 1e-12, "NS": 1e-9, "US": 1e-6, "MS": 1e-3, "S": 1, "M": 60, "H": 3600,
        "D": 86400, "Y": 3.15576e7}
HEADER = re.compile(r"#\s*\d+[A-Z]+\s*\(\s*([\d.]+(?:[Ee][+-]?\d+)?)\s*([A-Z]+)\s*\)")
n_h = n_hdr = 0
for (Z, A) in sorted(lib):
    with open(rdm / f"z{Z}.a{A}") as fh:
        h = HEADER.match(fh.readline())
    gs = [hl for E, flb, hl in ensdf.get((Z, A), []) if E == 0 and not flb and hl > 0]
    if not h or h.group(2) not in UNIT or not gs:
        continue
    n_hdr += 1
    t = float(h.group(1)) * UNIT[h.group(2)]
    if abs(t / gs[0] - 1) > tol_t:
        n_h += 1
        print(f"   {name(Z, A):16s} header {fmt_t(t):11s} ENSDFSTATE {fmt_t(gs[0]):11s} "
              f"ratio {t / gs[0]:.3g}   (z{Z}.a{A})")
print(f"   total: {n_h} of {n_hdr} headers compared")

sys.exit(1 if n_q or n_t or n_h else 0)
