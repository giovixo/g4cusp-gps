#!/usr/bin/env python3
"""
Report F: long-lived levels that Geant4 never decays correctly.

A level is decayed correctly only if
  (a) G4IonTable can create it, and
  (b) G4VRadioactiveDecay::LoadDecayTable finds decay data for it.
This script lists the levels for which (a) or (b) fails:

  1. levels without decay data: ENSDFSTATE levels with 1 s <= T1/2 < 1e15 s
     and no matching level in the decay library (matching rules of
     LoadDecayTable; excited levels without a floating base are excluded,
     because Geant4 de-excites them by IT anyway);
  2. shadowed E = 0 levels: at E = 0, G4IonTable::GetIon/FindIon and
     G4NuclideTable::GetIsotope ignore the floating-level base, so any request
     returns the first E = 0 level of the nuclide in ENSDFSTATE.dat, and the
     other E = 0 levels can never be created;
  3. decay records whose daughter is a shadowed level, i.e. decays whose
     daughter Geant4 replaces with another level, and what the replacement
     does: no decay data (radiation lost), other decay data, or the same data
     with another lifetime;
  4. what Geant4 does with each level of lists 1 and 2:
     A. created, no decay data: killed at its decay time without radiation
        (sub-case: an excited level, de-excited by IT instead of its own decay);
     B. not created, replaced by a level with decay data: wrong lifetime and
        radiation;
     C. not created, replaced by a level without decay data: no radiation.

Optionally it writes an rdecay01 macro that creates and decays every level of
list 1, and classifies the rdecay01 output of that macro.

Standard library only.
Usage:
    python3 check_undecayed_levels.py [RadioactiveDecay dir] [ENSDFSTATE.dat]
                                      [--min-parent=SECONDS] [--macro=FILE] [--log=FILE]
    --min-parent=S  list 3: only parents with T1/2 >= S (default 1 s)
    --macro=FILE    write an rdecay01 macro (200 single decays per level of list 1)
    --log=FILE      classify the output of rdecay01 run with that macro: each
                    level is "silent" (decays without any secondary), "not
                    created" (Geant4 gives back another level), or "decays"
                    (emits secondaries)
Defaults: $GEANT4_DATA_DIR/RadioactiveDecay6.1.2 and
          $GEANT4_DATA_DIR/G4ENSDFSTATE3.0/ENSDFSTATE.dat
"""
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

TOL = 1e-2          # keV, level-energy tolerance of LoadDecayTable (10 eV)
E0_TOL = 1e-3       # keV, G4NuclideTable level tolerance (1 eV)
TAU_MIN_NS = 1.0    # G4NuclideTable default mean-life threshold (1 ns)
LN2 = math.log(2)

data = Path(os.environ.get("GEANT4_DATA_DIR", "."))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
opts = dict((a[2:].split("=", 1) + [""])[:2] for a in sys.argv[1:] if a.startswith("--"))
rdm = Path(args[0]) if len(args) > 0 else data / "RadioactiveDecay6.1.2"
ensdf_file = Path(args[1]) if len(args) > 1 else data / "G4ENSDFSTATE3.0" / "ENSDFSTATE.dat"
min_parent = float(opts.get("min-parent") or 1.0)

SYM = ("n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge "
       "As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm "
       "Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U "
       "Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


def name(Z, A, E=0.0, flb=""):
    return f"{SYM[Z]}{A}" + (f"[{E:g}{flb}]" if E > 0 or flb else "")


def fmt_t(t):
    if t is None or t < 0:
        return "stable"
    for unit, sec in (("y", 3.15576e7), ("d", 86400), ("h", 3600), ("min", 60), ("s", 1)):
        if t >= sec:
            return f"{t / sec:.4g} {unit}"
    return f"{t:.3g} s"


def flag(tok):
    return "" if tok in ("-", "") else tok.lstrip("+-").upper()


# ---- ENSDFSTATE --------------------------------------------------------------------------
# ensdf:  (Z, A) -> [(E keV, flb, T1/2 s or -1)], all levels, file order
# loaded: the same, restricted to the levels G4NuclideTable loads (fixed ground states,
#         and levels with mean life >= 1 ns)
ensdf, loaded = defaultdict(list), defaultdict(list)
for line in open(ensdf_file):
    t = line.split()
    if len(t) < 6:
        continue
    Z, A, E, fl, tau = int(t[0]), int(t[1]), float(t[2]), flag(t[3]), float(t[4])
    lv = (E, fl, tau * 1e-9 * LN2 if tau > 0 else -1.0)
    ensdf[(Z, A)].append(lv)
    if (E == 0 and fl == "") or tau >= TAU_MIN_NS:
        loaded[(Z, A)].append(lv)

# ---- RadioactiveDecay: (Z, A) -> [parent level] in file order ----------------------------
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
            levels.append({"E": float(t[1]), "flb": flag(t[2]), "T": float(t[3]),
                           "summary": {}, "detail": []})
        elif levels and len(t) == 3:
            levels[-1]["summary"][t[0]] = levels[-1]["summary"].get(t[0], 0.0) + float(t[2])
        elif levels and len(t) >= 5:
            levels[-1]["detail"].append((t[0], float(t[1]), flag(t[2]), float(t[3])))
    lib[(int(m.group(1)), int(m.group(2)))] = levels


def g4_match(Z, A, E, flb):
    """Level chosen by G4VRadioactiveDecay::LoadDecayTable for an ion (Z, A, E, flb), or None."""
    for lv in lib.get((Z, A), []):
        if abs(lv["E"] - E) < TOL and (not flb or lv["flb"] == flb):
            return lv
    return None


def exact_match(Z, A, E, flb):
    for lv in lib.get((Z, A), []):
        if abs(lv["E"] - E) < TOL and lv["flb"] == flb:
            return lv
    return None


def e0_levels(Z, A):
    return [lv for lv in loaded.get((Z, A), []) if abs(lv[0]) < E0_TOL]


print(f"RadioactiveDecay: {rdm}\nENSDFSTATE:       {ensdf_file}\n")

# ---- 1. levels without decay data ----------------------------------------------------------
print("1. ENSDFSTATE levels with 1 s <= T1/2 < 1e15 s and no decay data in the library")
missing = []
for (Z, A), levels in sorted(ensdf.items()):
    for E, flb, hl in levels:
        if 1.0 <= hl < 1e15 and g4_match(Z, A, E, flb) is None and not (E > 0 and not flb):
            missing.append((Z, A, E, flb, hl))
            where = "no file" if (Z, A) not in lib else "no level"
            print(f"   {name(Z, A, E, flb):16s} T1/2 = {fmt_t(hl):10s} ({where})")
print(f"   total: {len(missing)} (excited levels without flag are excluded: Geant4 applies IT to them)\n")

# ---- 2. shadowed E = 0 levels ---------------------------------------------------------------
multi = {k: e0_levels(*k) for k in loaded if len(e0_levels(*k)) > 1}
shadowed = [(k, lv) for k, lvs in multi.items() for lv in lvs[1:]]
print(f"2. Nuclides with more than one E = 0 level loaded by G4NuclideTable: {len(multi)}")
print(f"   Shadowed E = 0 levels (never created; any E = 0 request returns the first level): "
      f"{len(shadowed)}")
long = [(k, lv) for k, lv in shadowed if lv[2] >= 1.0]
n_data = sum(1 for (Z, A), lv in long if exact_match(Z, A, 0.0, lv[1]))
print(f"   of which with T1/2 >= 1 s: {len(long)}, {n_data} of them with their own decay data "
      f"in the library (never used)")
print(f"   First level is the fixed ground state for "
      f"{sum(1 for lvs in multi.values() if lvs[0][1] == '')} nuclides, a floating level for "
      f"{sum(1 for lvs in multi.values() if lvs[0][1] != '')}.\n")
print(f"   {'shadowed level':18s} {'T1/2':>11s}  data  {'Geant4 creates instead':22s} {'T1/2':>11s}  data")
for (Z, A), lv in sorted(long, key=lambda s: -s[1][2]):
    g = multi[(Z, A)][0]
    print(f"   {name(Z, A, 0, lv[1] or '-'):18s} {fmt_t(lv[2]):>11s}  "
          f"{'yes' if exact_match(Z, A, 0.0, lv[1]) else 'NO ':4s}  "
          f"{name(Z, A, 0, g[1] or '-'):22s} {fmt_t(g[2]):>11s}  "
          f"{'yes' if g4_match(Z, A, 0.0, g[1]) else 'NO'}")
print()

# ---- 3. decay records feeding a shadowed level ----------------------------------------------
DAUGHTER = {"BetaMinus": (1, 0), "BetaPlus": (-1, 0), "KshellEC": (-1, 0), "LshellEC": (-1, 0),
            "MshellEC": (-1, 0), "NshellEC": (-1, 0), "Alpha": (-2, -4), "Proton": (-1, -1),
            "Neutron": (0, -1), "IT": (0, 0)}
merged, info = defaultdict(float), {}
shadow_parent = {(Z, A, lv[1]) for (Z, A), lv in shadowed}
skipped_parents = set()
for (Z, A), levels in lib.items():
    hl_ensdf = {(round(E, 3), fl): T for E, fl, T in ensdf.get((Z, A), [])}
    for lv in levels:
        tP = hl_ensdf.get((round(lv["E"], 3), lv["flb"]), lv["T"])
        if tP is None or tP < min_parent:
            continue
        if abs(lv["E"]) <= E0_TOL and (Z, A, lv["flb"]) in shadow_parent:
            skipped_parents.add(name(Z, A, 0, lv["flb"]))   # parent never created: record unused
            continue
        psum = defaultdict(float)
        for mode, Ed, fd, pct in lv["detail"]:
            psum[mode] += pct
        for mode, Ed, fd, pct in lv["detail"]:
            if mode not in DAUGHTER or abs(Ed) > E0_TOL:
                continue
            Zd, Ad = Z + DAUGHTER[mode][0], A + DAUGHTER[mode][1]
            lvs = e0_levels(Zd, Ad)
            if len(lvs) < 2 or lvs[0][1] == fd:
                continue
            req = next((x for x in lvs if x[1] == fd), None)
            if req is None:
                continue
            br = lv["summary"].get(mode, 0.0) * pct / psum[mode] if psum[mode] > 0 else 0.0
            if br <= 0:
                continue
            g = lvs[0]
            d_req, d_g4 = exact_match(Zd, Ad, 0.0, fd), g4_match(Zd, Ad, 0.0, g[1])
            if d_g4 is None:
                cls = "no radiation" if d_req else "no data for either level"
            elif d_req is not None and d_g4 is d_req:
                cls = "same data, other lifetime"
            else:
                cls = "other lifetime and decay data"
            key = (Z, A, lv["E"], lv["flb"], Zd, Ad, fd)
            merged[key] += br              # EC sub-shells and beta+ merged
            info[key] = (tP, req[2], g, cls)
classes = defaultdict(int)
for k in merged:
    classes[info[k][3]] += 1
print(f"3. Decays (parent T1/2 >= {min_parent:g} s) whose daughter is a shadowed E = 0 level: "
      f"{len(merged)} parent -> daughter pairs")
for c, n in sorted(classes.items(), key=lambda x: -x[1]):
    print(f"   {c}: {n}")
if skipped_parents:
    print(f"   (parents that are themselves shadowed E = 0 levels are skipped: their records are "
          f"never used; {len(skipped_parents)} with T1/2 >= {min_parent:g} s)")
print(f"\n   {'parent':14s} {'T1/2':>10s}  {'requested daughter':18s} {'T1/2':>10s}  "
      f"{'Geant4 daughter':16s} {'T1/2':>10s}  {'BR':>7s}  consequence")
for k in sorted(merged, key=lambda k: (-merged[k] * min(info[k][0], 3.15576e9), k)):
    Z, A, E, fl, Zd, Ad, fd = k
    tP, tR, g, cls = info[k]
    print(f"   {name(Z, A, E, fl):14s} {fmt_t(tP):>10s}  {name(Zd, Ad, 0, fd or '-'):18s} "
          f"{fmt_t(tR):>10s}  {name(Zd, Ad, 0, g[1] or '-'):16s} {fmt_t(g[2]):>10s}  "
          f"{100 * merged[k]:6.2f}%  {cls}")

# ---- 4. what Geant4 does with the levels of lists 1 and 2 ----------------------------------
# Predicted from the data and the rules of G4IonTable and LoadDecayTable; the rdecay01 run
# (--macro, --log) checks the levels of list 1.
shadow_long = {(Z, A, lv[1]): (Z, A, lv) for (Z, A), lv in long}
CAT_A = "A. created; no decay data"
CAT_B = "B. not created; replaced by a level with decay data: wrong lifetime and radiation"
CAT_C = "C. not created; replaced by a level without decay data: no radiation"
outcome = {CAT_A: [], CAT_B: [], CAT_C: []}
it_case = []        # sub-case of A: excited level de-excited by IT instead of its own decay
union = {(Z, A, E, flb): hl for Z, A, E, flb, hl in missing}
for (Z, A), lv in long:
    union.setdefault((Z, A, 0.0, lv[1]), lv[2])
for (Z, A, E, flb), hl in sorted(union.items()):
    label = f"{name(Z, A, E, flb)} ({fmt_t(hl)})"
    if E == 0 and (Z, A, flb) in shadow_long:
        g = multi[(Z, A)][0]
        own = "own data unused" if exact_match(Z, A, 0.0, flb) else "no own data"
        sub = f"{label} -> {name(Z, A, 0, g[1] or '-')} ({fmt_t(g[2])}), {own}"
        outcome[CAT_C if g4_match(Z, A, 0.0, g[1]) is None else CAT_B].append(sub)
    else:
        outcome[CAT_A].append(label)
        if E > 0:
            it_case.append(label)
n_a, n_b, n_c = (len(outcome[k]) for k in (CAT_A, CAT_B, CAT_C))
print(f"\n4. What Geant4 does with the {len(union)} levels of lists 1 and 2 "
      f"({len(missing)} + {len(long)} - {len(missing) + len(long) - len(union)} in both)")
for k in outcome:
    print(f"   {k}: {len(outcome[k])}")
    if k == CAT_A:
        print(f"      killed at their decay time without radiation: {n_a - len(it_case)}")
        print(f"      de-excited by IT instead of their own decay: {len(it_case)}")
print(f"   no radiation at all (A except IT, and C): {n_a - len(it_case) + n_c}")
print(f"   radiation of the wrong level or decay mode (IT and B): {len(it_case) + n_b}")
for k in outcome:
    if k == CAT_A and len(outcome[k]) > 20:
        if it_case:
            print(f"\n   {k}, de-excited by IT instead of their own decay (sub-case of A):")
            for s in it_case:
                print(f"      {s}")
    elif len(outcome[k]) <= 20:
        print(f"\n   {k}:")
        for s in outcome[k]:
            print(f"      {s}")

# ---- rdecay01 macro and classification (list 1) --------------------------------------------
if opts.get("macro"):
    with open(opts["macro"], "w") as f:
        f.write("# Generated by check_undecayed_levels.py: one run of 200 single decays per level\n"
                "# of list 1 (levels without decay data).\n"
                "# Run with the unmodified example examples/extended/radioactivedecay/rdecay01\n"
                "# and classify the output with check_undecayed_levels.py --log=<output>.\n"
                "/control/verbose 0\n/run/verbose 0\n"
                "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year\n"
                "/run/initialize\n/gun/particle ion\n/rdecay01/fullChain false\n/run/printProgress 0\n")
        for Z, A, E, flb, hl in missing:
            f.write(f"/control/echo CASE {Z} {A} {E:g} {flb or '-'}\n"
                    f"/gun/ion {Z} {A} 0 {E:g} {flb}\n/run/beamOn 200\n")
    print(f"\nrdecay01 macro written to {opts['macro']}")

if opts.get("log"):
    # Split the output at the "CASE" echo lines; in each block read the primary name
    # ("The run was N <name> of") and the "Nb of generated particles" table.
    cases, cur = [], None
    for line in open(opts["log"], errors="replace"):
        if line.startswith("CASE "):
            cur = {"req": line.split()[1:], "primary": None, "rows": [], "table": False}
            cases.append(cur)
        elif cur is not None:
            m = re.match(r"\s*The run was \d+ (\S+) of", line)
            if m:
                cur["primary"] = m.group(1)
            elif "Nb of generated particles" in line:
                cur["table"] = True
            elif cur["table"]:
                m = re.match(r"\s*(\S+):\s+\d+\s+Emean", line)
                if m:
                    cur["rows"].append(m.group(1))
                elif cur["rows"] and not line.strip():
                    cur["table"] = False
    shadow_set = {(Z, A, lv[1]) for (Z, A), lv in shadowed}
    out = {"silent": [], "not created": [], "decays": []}
    for c in cases:
        Z, A, E, flb = int(c["req"][0]), int(c["req"][1]), float(c["req"][2]), c["req"][3].strip("-")
        m = re.fullmatch(r"([A-Za-z]+\d+)(?:\[([\d.]+)([A-Z]?)\])?", c["primary"] or "")
        got_flb = (m.group(3) or "") if m else None
        got_E = float(m.group(2)) if m and m.group(2) else 0.0
        hl = next((h for z, a, e, f, h in missing if (z, a, e, f) == (Z, A, E, flb)), -1)
        label = f"{name(Z, A, E, flb):16s} T1/2 = {fmt_t(hl):10s} Geant4 primary: {c['primary']}"
        if m is None or got_flb != flb or abs(got_E - E) > TOL:
            sh = "  (shadowed, list 2)" if E == 0 and (Z, A, flb) in shadow_set else ""
            out["not created"].append(label + sh)
        elif len(c["rows"]) <= 1:
            out["silent"].append(label)
        else:
            out["decays"].append(label + "  secondaries: "
                                 + ", ".join([r for r in c["rows"] if r != c["primary"]][:6]))
    print(f"\nrdecay01 results from {opts['log']} ({len(cases)} levels)")
    for k in ("silent", "not created", "decays"):
        print(f"\n  {k}: {len(out[k])}")
        for s in out[k]:
            print("   ", s)
