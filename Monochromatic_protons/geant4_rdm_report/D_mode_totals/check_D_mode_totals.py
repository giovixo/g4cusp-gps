#!/usr/bin/env python3
"""
Report D: candidate wrong decay-mode totals in RadioactiveDecay6.1.2.

G4VRadioactiveDecay::LoadDecayTable rescales the detail lines of each mode so
that they add up to the mode's summary fraction: only the summary lines decide
how often each mode occurs.  This script lists levels (1 s <= T1/2 < 1e15 s)
where the summary looks inconsistent with the detail lines:
  - the summary gives > 50 % of the decays to a mode whose detail lines add up
    to < 50 %, or
  - the summary gives ~0 to a mode whose detail lines add up to >= 1 %.
This is a heuristic: each candidate has to be checked against ENSDF.

Standard library only.
Usage:
    python3 check_D_mode_totals.py [RadioactiveDecay dir] [ENSDFSTATE.dat]
Defaults: $GEANT4_DATA_DIR/RadioactiveDecay6.1.2 and
          $GEANT4_DATA_DIR/G4ENSDFSTATE3.0/ENSDFSTATE.dat
"""
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

TOL = 1e-2   # keV, level-energy tolerance

data = Path(os.environ.get("GEANT4_DATA_DIR", "."))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
rdm = Path(args[0]) if len(args) > 0 else data / "RadioactiveDecay6.1.2"
ensdf_file = Path(args[1]) if len(args) > 1 else data / "G4ENSDFSTATE3.0" / "ENSDFSTATE.dat"

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
ensdf = defaultdict(list)
for line in open(ensdf_file):
    t = line.split()
    if len(t) >= 5:
        tau_ns = float(t[4])
        ensdf[(int(t[0]), int(t[1]))].append(
            (float(t[2]), t[3].lstrip("+-").upper(), tau_ns * 1e-9 * math.log(2) if tau_ns > 0 else -1))

# ---- RadioactiveDecay: (Z, A) -> [level dict] in file order --------------------
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


def g4_match(Z, A, E, flb):
    """Level chosen by G4VRadioactiveDecay::LoadDecayTable, or None."""
    for lv in lib.get((Z, A), []):
        if abs(lv["E"] - E) < TOL and (not flb or lv["flb"] == flb):
            return lv
    return None


print(f"RadioactiveDecay: {rdm}\nENSDFSTATE:       {ensdf_file}\n")

NOT_IMPL = {"BDProton", "BDNeutron", "Beta2Minus", "Beta2Plus", "Proton2", "Neutron2"}
print("Candidate wrong decay-mode totals (summary > 0.5 with detail lines < 50 %, "
      "or summary ~0 with detail lines >= 1 %)")
n = 0
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        hl = [h for e, f, h in ensdf.get((Z, A), []) if abs(e - lv["E"]) < TOL and f == lv["flb"]]
        if not hl or not 1.0 <= hl[0] < 1e15:
            continue
        dsum = defaultdict(float)
        for d in lv["detail"]:
            dsum[d[0]] += d[2]
        bad = [m for m, v in lv["summary"].items()
               if m != "IT" and m not in NOT_IMPL and v > 0.5 and dsum.get(m, 0) < 50]
        bad += [m for m, v in dsum.items()          # mode with intensity but no total
                if m not in NOT_IMPL and v >= 1 and lv["summary"].get(m, 0) < 1e-6]
        if bad:
            n += 1
            summ = ", ".join(f"{m} {v:.3g}" for m, v in lv["summary"].items() if v > 1e-6)
            det = ", ".join(f"{m} {v:.3g}%" for m, v in dsum.items() if v > 0)
            print(f"   {name(Z, A, lv['E'], lv['flb']):16s} T1/2 = {fmt_t(hl[0]):10s} "
                  f"summary: {summ}; detail sums: {det}")
print(f"   total: {n}")
