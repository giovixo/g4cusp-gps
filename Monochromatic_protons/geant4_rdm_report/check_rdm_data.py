#!/usr/bin/env python3
"""
Consistency checks of the Geant4 RadioactiveDecay data library against
G4ENSDFSTATE.  Standard library only.

Usage:
    python3 check_rdm_data.py [RadioactiveDecay dir] [ENSDFSTATE.dat]
Defaults: $GEANT4_DATA_DIR/RadioactiveDecay6.1.2 and
          $GEANT4_DATA_DIR/G4ENSDFSTATE3.0/ENSDFSTATE.dat

Checks
  A. ground states with an IT branch (decay onto themselves)
  B. decay branches (non-IT) with Q <= 0 in the detail lines
  C. levels with 1 s <= T1/2 < 1e15 s in ENSDFSTATE but no matching level in
     the decay library (matching rules of G4VRadioactiveDecay::LoadDecayTable)
  D. ground states without floating flag in ENSDFSTATE that LoadDecayTable
     matches to a floating level of the decay library
  E. 'P'-line half-lives differing by > 1 % from ENSDFSTATE (informational:
     LoadDecayTable ignores them)
"""
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

TOL = 1e-2   # keV

data = Path(os.environ.get("GEANT4_DATA_DIR", "."))
rdm = Path(sys.argv[1]) if len(sys.argv) > 1 else data / "RadioactiveDecay6.1.2"
ensdf_file = Path(sys.argv[2]) if len(sys.argv) > 2 else data / "G4ENSDFSTATE3.0" / "ENSDFSTATE.dat"

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

# ---- A ----------------------------------------------------------------------
print("A. Ground states (no floating flag) with an IT branch onto themselves")
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        if lv["E"] == 0 and not lv["flb"] and lv["summary"].get("IT", 0) > 0:
            others = ", ".join(f"{m} {v:g}" for m, v in lv["summary"].items() if m != "IT")
            print(f"   {name(Z, A):8s} IT {lv['summary']['IT']:<8g} (other modes: {others})")

# ---- B ----------------------------------------------------------------------
print("\nB. Levels with non-IT decay branches having Q <= 0 (intensity > 0)")
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        bad = [d for d in lv["detail"] if d[0] != "IT" and d[2] > 0 and d[3] <= 0]
        if bad:
            frac = sum(d[2] for d in bad)
            modes = sorted({d[0] for d in bad})
            qs = ", ".join(f"{d[3]:g}" for d in bad[:3])
            print(f"   {name(Z, A, lv['E'], lv['flb']):16s} {len(bad):3d} branches, {frac:8.3f} % "
                  f"of decays, {'/'.join(modes)}; Q = {qs}{', ...' if len(bad) > 3 else ''}")

# ---- C ----------------------------------------------------------------------
print("\nC. ENSDFSTATE levels with 1 s <= T1/2 < 1e15 s and no decay data in the library")
n_c = 0
for (Z, A), levels in sorted(ensdf.items()):
    for E, flb, hl in levels:
        if 1.0 <= hl < 1e15 and g4_match(Z, A, E, flb) is None and not (E > 0 and not flb):
            n_c += 1
            where = "no file" if (Z, A) not in lib else "no level"
            print(f"   {name(Z, A, E, flb):16s} T1/2 = {fmt_t(hl):10s} ({where})")
print(f"   total: {n_c} (excited levels without flag are excluded: Geant4 applies IT to them)")

# ---- D ----------------------------------------------------------------------
print("\nD. ENSDFSTATE levels without floating flag that LoadDecayTable matches to a floating level")
for (Z, A), levels in sorted(ensdf.items()):
    for E, flb, hl in levels:
        if flb or hl < 1.0:
            continue
        exact = any(abs(lv["E"] - E) < TOL and not lv["flb"] for lv in lib.get((Z, A), []))
        lv = g4_match(Z, A, E, flb)
        if lv is not None and not exact:
            ens_float = [h for e, f, h in levels if abs(e - E) < TOL and f == lv["flb"]]
            modes = ", ".join(f"{m} {v:g}" for m, v in lv["summary"].items())
            print(f"   {name(Z, A, E):14s} T1/2 = {fmt_t(hl):9s} -> library level "
                  f"{name(Z, A, lv['E'], lv['flb'])} (T1/2 = "
                  f"{fmt_t(ens_float[0]) if ens_float else '?'}): {modes}")

# ---- E ----------------------------------------------------------------------
print("\nE. 'P'-line half-lives differing from ENSDFSTATE by more than 1 %")
n_e = 0
for (Z, A), levels in sorted(lib.items()):
    for lv in levels:
        if lv["T"] <= 0:
            continue
        for E, flb, hl in ensdf.get((Z, A), []):
            if abs(E - lv["E"]) < TOL and flb == lv["flb"]:
                if hl > 0 and abs(lv["T"] / hl - 1) > 0.01:
                    n_e += 1
                    print(f"   {name(Z, A, lv['E'], lv['flb']):16s} P line {fmt_t(lv['T']):11s} "
                          f"ENSDFSTATE {fmt_t(hl)}")
                break
print(f"   total: {n_e}")
