#!/usr/bin/env python3
"""
Report E: RadioactiveDecay6.1.2 'P'-line half-lives and README_RDM format.

  1. 'P'-line half-lives that differ by more than 1 % from ENSDFSTATE
     (G4VRadioactiveDecay::LoadDecayTable ignores them).
  2. Actual format of the decay files, to compare with README_RDM:
     column counts of the summary and detail lines, the line-length rule that
     LoadDecayTable uses to tell them apart (after removing trailing blanks:
     summary < 72 characters, detail without beta type < 84, with beta type
     >= 84; lines are read into a 120-character buffer), and the convention
     of the detail percentages (relative to the mode total, as README_RDM
     says, or relative to all decays).

Standard library only.
Usage:
    python3 check_E_halflives_readme.py [RadioactiveDecay dir] [ENSDFSTATE.dat]
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

print("1. 'P'-line half-lives differing from ENSDFSTATE by more than 1 %")
n = 0
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        if lv["T"] <= 0:
            continue
        for E, flb, hl in ensdf.get((Z, A), []):
            if abs(E - lv["E"]) < TOL and flb == lv["flb"]:
                if hl > 0 and abs(lv["T"] / hl - 1) > 0.01:
                    n += 1
                    print(f"   {name(Z, A, lv['E'], lv['flb']):16s} P line {fmt_t(lv['T']):11s} "
                          f"ENSDFSTATE {fmt_t(hl)}")
                break
print(f"   total: {n}")

print("\n2. File format")
cols = defaultdict(int)
rule_bad, longest = [], 0
for f in sorted(rdm.glob("z*.a*")):
    for line in f.read_text().splitlines():
        longest = max(longest, len(line))
        s = line.rstrip()
        t = s.split()
        if not t or s.startswith("#") or s.startswith("P"):
            continue
        cols[len(t)] += 1
        ok = (len(t) == 3 and len(s) < 72) or (len(t) == 5 and 72 <= len(s) < 84) \
            or (len(t) == 6 and len(s) >= 84)
        if not ok:
            rule_bad.append(f"{f.name}: {s}")
print("   decay-mode lines by number of columns: "
      + ", ".join(f"{k} columns: {v}" for k, v in sorted(cols.items())))
print(f"   lines breaking the length rule of LoadDecayTable: {len(rule_bad)}")
for s in rule_bad[:20]:
    print("     ", s)
print(f"   longest line: {longest} characters (buffer: 120)")

n_det = n_all = 0
for levels in lib.values():
    for lv in levels:
        if not lv["detail"]:
            continue
        n_det += 1
        dsum = defaultdict(float)
        for d in lv["detail"]:
            dsum[d[0]] += d[2]
        if any(abs(v - 100) > 1 for v in dsum.values()):
            n_all += 1
print(f"   levels with detail lines: {n_det}; levels where the detail percentages of some mode "
      f"do not add up to 100 (+-1): {n_all}")
