"""
decay_chain_builder.py
======================
Build and serialise radioactive decay chains from the Geant4 data libraries.

Usage (CLI)
-----------
    python decay_chain_builder.py Os172 --outdir DecayChains/
    python decay_chain_builder.py Pm139 --prune 0.001 --min-halflife 1e-6
    python decay_chain_builder.py Re182 --db /path/to/RadioactiveDecay6.1.2

Usage (library)
---------------
    from decay_chain_builder import DecayDatabase, build_chain, write_chain

    db    = DecayDatabase()                 # locate the Geant4 data automatically
    chain = build_chain("Os172", db)
    write_chain("Os172", chain, "DecayChains/")

Requirements: Python 3.9+, no third-party libraries.

Data sources
------------
Decay modes and branching ratios come from the RadioactiveDecay library
(files z<Z>.a<A>).  Half-lives come from ENSDFSTATE.dat (G4ENSDFSTATE library),
exactly as in Geant4: the README of RadioactiveDecay states that the half-life
on its 'P' lines is ignored, and some of those values are wrong (e.g. ground
state and isomer swapped for Sc44, Eu152, Cu68; Se79 off by ten orders of
magnitude).  The 'P'-line value is used only for a level missing from
ENSDFSTATE.

Database file format (README_RDM)
---------------------------------
    P  <E keV>  <flag>  <T1/2>              parent level; flag '-' or '+X', '+Y', ...
       <mode>  0  <fraction>                summary: total branching of the mode
       <mode>  <Ed keV>  <flag>  <%>  <Q> [<forbiddenness>]
                                            detail: daughter level Ed; % of all decays
Each mode's total comes from its summary line.  It is shared among the
daughter levels in proportion to the detail lines (which do not always add up
to the summary exactly); a mode without detail lines feeds the daughter
ground state.

Levels are identified by (Z, A, E, floating-level letter) and matched to the
library exactly as Geant4 does (see DecayDatabase.get_level), so that the
chains describe what the Geant4 decay simulation will do.  A warning is
printed for every level whose data had to be borrowed or assumed.  An excited
level absent from both libraries de-excites promptly to the ground state.

Output format (one file per root isotope, e.g. DecayChains/Os172.dat)
----------------------------------------------------------------------
    # Os172 Decay Chain
    # params: format=4 prune=0.001 min_halflife=1e-06 db=RadioactiveDecay6.1.2 ensdf=G4ENSDFSTATE3.0
    # HalfLife[s] (Z_daughter, A_daughter, E_daughter[, flb])
    # Chain n. 0 Branching ratio: 0.986200
    1.920000e+01 (75, 172, 0.0)
    ...
Each data line is one decay STEP: the half-life of the parent nucleus, then
the daughter it produces.  The floating-level letter is written only when set.

Design notes
------------
* Database files are parsed lazily and cached.
* Loop detection uses the set of nuclides already on the current path.
* Nuclides with half-life below min_halflife are not recorded: a step goes
  from a long-lived parent directly to the first long-lived (or stable)
  nuclides reached through such prompt decays, with the product of the
  branching ratios.  The daughter of step i is always the parent of step i+1.
* Transitions of a level onto itself (spurious IT branches of some ground
  states in RadioactiveDecay6.1.2) are dropped and the branching renormalised.
* Branches whose cumulative branching ratio falls below prune_threshold are
  dropped (default 0.001, as in Campana et al. 2026); the remaining chains are
  renormalised to a total of 1.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path
from typing import NamedTuple

from nuclides import FLOAT_LEVEL_CHARS, SYMBOLS, format_name, parse_name

# Backward-compatible aliases (used by compute_activities.py)
_INV_PERIODIC_TABLE = SYMBOLS

_E_TOL = 1e-2          # keV, tolerance when matching level energies
CHAIN_FORMAT = 4       # bump when the chain-building algorithm changes (forces rebuilds)
_LN2 = math.log(2)


# ---------------------------------------------------------------------------
# Nuclide identity type
# ---------------------------------------------------------------------------

class Nuclide(NamedTuple):
    Z: int
    A: int
    E: float        # excitation energy [keV]
    flb: str = ""   # floating level letter ('' = fixed level)

    def __str__(self) -> str:
        return format_name(self.Z, self.A, self.E, self.flb)


def _flag_to_flb(flag: str) -> str:
    """'-' -> '', '+X' -> 'X'."""
    flb = flag.lstrip("+-").upper()
    if flb and flb not in FLOAT_LEVEL_CHARS:
        raise ValueError(f"unknown floating-level flag {flag!r}")
    return flb


# ---------------------------------------------------------------------------
# DecayStep: one edge in the decay tree
# ---------------------------------------------------------------------------

class DecayStep(NamedTuple):
    parent_halflife: float      # half-life of the decaying nucleus [s]
    daughter: Nuclide           # identity of the daughter nucleus
    branch_ratio: float         # branching ratio for THIS edge (0–1)


# A complete chain is a list of DecayStep from root to terminal.
Chain = list[DecayStep]


class LevelData(NamedTuple):
    """Decay data for one level of a nuclide."""
    Z: int
    A: int
    excitation: float               # keV
    flb: str                        # floating level letter
    half_life: float                # seconds (< 0 = stable)
    daughters: list[tuple[Nuclide, float]]  # (daughter, branching ratio), sum = 1
    synthetic: bool = False         # True: not in the decay library (IT assumed)


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

_DELTA: dict[str, tuple[int, int]] = {
    "Alpha":      (-2, -4),
    "BetaPlus":   (-1,  0),
    "BetaMinus":  (+1,  0),
    "IT":         ( 0,  0),
    "Neutron":    ( 0, -1),
    "Proton":     (-1, -1),
    "Beta2Minus": (+2,  0),
    "Beta2Plus":  (-2,  0),
    "Triton":     (-1, -3),
}


def _mode_delta(mode: str) -> tuple[int, int] | None:
    if mode in _DELTA:
        return _DELTA[mode]
    if mode.endswith("EC"):             # KshellEC, LshellEC, MshellEC, NshellEC
        return (-1, 0)
    return None


def _parse_decay_file(path: Path) -> list[LevelData]:
    """
    Parse one RadioactiveDecay file into LevelData (half-lives from the 'P'
    lines; DecayDatabase replaces them with the ENSDFSTATE values).
    """
    m = re.fullmatch(r'z(\d+)\.a(\d+)', path.name)
    if not m:
        raise ValueError(f"Unexpected filename format: {path.name}")
    Z, A = int(m.group(1)), int(m.group(2))

    blocks: list[dict] = []
    for raw in path.read_text().splitlines():
        tok = raw.split()
        if not tok or tok[0].startswith("#"):
            continue
        if tok[0] == "P":
            blocks.append({"E": float(tok[1]), "flb": _flag_to_flb(tok[2]),
                           "T": float(tok[3]), "summary": {}, "detail": {}})
            continue
        if not blocks:
            continue
        b, mode = blocks[-1], tok[0]
        if len(tok) == 3:                                   # summary line
            b["summary"][mode] = b["summary"].get(mode, 0.0) + float(tok[2])
        elif len(tok) in (5, 6):                            # detail line
            key = (float(tok[1]), _flag_to_flb(tok[2]))
            d = b["detail"].setdefault(mode, {})
            d[key] = d.get(key, 0.0) + float(tok[3])
        else:
            print(f"  [WARNING] unexpected line in {path.name}: {raw.strip()!r}",
                  file=sys.stderr)

    levels = []
    for b in blocks:
        acc: dict[Nuclide, float] = {}
        modes = set(b["summary"]) | set(b["detail"])
        for mode in modes:
            if mode == "SpFission":
                continue                    # fission fragments: not followed
            delta = _mode_delta(mode)
            if delta is None:
                print(f"  [WARNING] unknown decay mode '{mode}' in {path.name}, skipped",
                      file=sys.stderr)
                continue
            dz, da = delta
            detail = b["detail"].get(mode, {})
            det_sum = sum(detail.values())
            total = b["summary"].get(mode, det_sum / 100.0)
            if total <= 0:
                continue
            if det_sum > 0:
                for (ed, flb), pct in detail.items():
                    d = Nuclide(Z + dz, A + da, ed, flb)
                    acc[d] = acc.get(d, 0.0) + total * pct / det_sum
            else:
                d = Nuclide(Z + dz, A + da, 0.0, "")
                acc[d] = acc.get(d, 0.0) + total
        # Drop transitions of a level onto itself: ten ground states of
        # RadioactiveDecay6.1.2 (e.g. Y92, Tb154, Ir186) carry a spurious IT
        # branch, which cannot occur from a ground state.
        parent = Nuclide(Z, A, b["E"], b["flb"])
        acc = {d: br for d, br in acc.items()
               if not (d.Z == Z and d.A == A and abs(d.E - parent.E) < _E_TOL
                       and d.flb == parent.flb)}
        norm = sum(acc.values())
        daughters = [(d, br / norm) for d, br in acc.items()] if norm > 0 else []
        levels.append(LevelData(Z, A, b["E"], b["flb"], b["T"], daughters))
    return levels


def _parse_ensdf(path: Path) -> dict[tuple[int, int], list[tuple[float, str, float]]]:
    """ENSDFSTATE.dat -> {(Z, A): [(E keV, flb, half-life s or -1 if stable)]}."""
    table: dict[tuple[int, int], list[tuple[float, str, float]]] = {}
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) < 5:
                continue
            tau_ns = float(tok[4])
            hl = tau_ns * 1e-9 * _LN2 if tau_ns > 0 else -1.0
            table.setdefault((int(tok[0]), int(tok[1])), []).append(
                (float(tok[2]), _flag_to_flb(tok[3]), hl))
    return table


# ---------------------------------------------------------------------------
# Locating the Geant4 data
# ---------------------------------------------------------------------------

def _latest(parent: Path, prefix: str) -> Path | None:
    cands = sorted(p for p in parent.glob(prefix + "*") if p.is_dir())
    return cands[-1] if cands else None


def _resolve_db(db_arg: str | None) -> Path:
    """RadioactiveDecay directory: --db, $G4RADIOACTIVEDATA, $GEANT4_DATA_DIR, ./"""
    if db_arg:
        return Path(db_arg)
    if os.environ.get("G4RADIOACTIVEDATA"):
        return Path(os.environ["G4RADIOACTIVEDATA"])
    for parent in (os.environ.get("GEANT4_DATA_DIR"), "."):
        if parent and (p := _latest(Path(parent), "RadioactiveDecay")):
            return p
    raise FileNotFoundError(
        "Cannot locate the Geant4 RadioactiveDecay library: use --db, or set "
        "$G4RADIOACTIVEDATA or $GEANT4_DATA_DIR (source geant4.sh)."
    )


def _resolve_ensdf(ensdf_arg: str | None) -> Path:
    """ENSDFSTATE.dat: --ensdf, $G4ENSDFSTATEDATA, $GEANT4_DATA_DIR, ./"""
    if ensdf_arg:
        p = Path(ensdf_arg)
        return p / "ENSDFSTATE.dat" if p.is_dir() else p
    if os.environ.get("G4ENSDFSTATEDATA"):
        return Path(os.environ["G4ENSDFSTATEDATA"]) / "ENSDFSTATE.dat"
    for parent in (os.environ.get("GEANT4_DATA_DIR"), "."):
        if parent and (p := _latest(Path(parent), "G4ENSDFSTATE")):
            return p / "ENSDFSTATE.dat"
    raise FileNotFoundError(
        "Cannot locate the Geant4 ENSDFSTATE library: use --ensdf, or set "
        "$G4ENSDFSTATEDATA or $GEANT4_DATA_DIR (source geant4.sh)."
    )


# ---------------------------------------------------------------------------
# DecayDatabase
# ---------------------------------------------------------------------------

class DecayDatabase:
    """
    Decay data from the Geant4 RadioactiveDecay library, with half-lives from
    ENSDFSTATE.  Files are parsed on demand and cached.
    """

    def __init__(self, db_path: str | Path | None = None,
                 ensdf_path: str | Path | None = None, verbose: bool = True):
        self._db_path = _resolve_db(str(db_path) if db_path else None)
        if not self._db_path.is_dir():
            raise FileNotFoundError(f"Database directory not found: {self._db_path}")
        self._ensdf_path = _resolve_ensdf(str(ensdf_path) if ensdf_path else None)
        if not self._ensdf_path.is_file():
            raise FileNotFoundError(f"ENSDFSTATE file not found: {self._ensdf_path}")

        self._index: dict[tuple[int, int], Path] = {}
        for f in self._db_path.iterdir():
            m = re.fullmatch(r'z(\d+)\.a(\d+)', f.name)
            if m:
                self._index[(int(m.group(1)), int(m.group(2)))] = f
        if not self._index:
            raise FileNotFoundError(f"No z*.a* files in {self._db_path}")
        self._ensdf = _parse_ensdf(self._ensdf_path)
        self._cache: dict[tuple[int, int], list[LevelData]] = {}
        self.verbose = verbose
        self.n_hl_corrected = 0
        self._warned: set[Nuclide] = set()

    # -- identification of the data used (written in the chain files) --
    @property
    def tag(self) -> str:
        return f"db={self._db_path.name} ensdf={self._ensdf_path.parent.name}"

    @property
    def size(self) -> int:
        return len(self._index)

    def has_file(self, Z: int, A: int) -> bool:
        """True if the decay library has a file for (Z, A)."""
        return (Z, A) in self._index

    def ensdf_halflife(self, nuc: Nuclide) -> float | None:
        """ENSDFSTATE half-life [s] of a level (-1 = stable), or None if absent."""
        for e, flb, hl in self._ensdf.get((nuc.Z, nuc.A), []):
            if abs(e - nuc.E) < _E_TOL and flb == nuc.flb:
                return hl
        return None

    def get_levels(self, Z: int, A: int) -> list[LevelData]:
        """All levels of (Z, A) in the decay library, half-lives from ENSDFSTATE."""
        key = (Z, A)
        if key not in self._cache:
            path = self._index.get(key)
            levels = _parse_decay_file(path) if path else []
            fixed = []
            for lv in levels:
                hl = self.ensdf_halflife(Nuclide(Z, A, lv.excitation, lv.flb))
                if hl is not None and hl != lv.half_life:
                    if lv.half_life > 0 and abs(hl / lv.half_life - 1) > 0.01:
                        self.n_hl_corrected += 1
                    lv = lv._replace(half_life=hl)
                fixed.append(lv)
            self._cache[key] = fixed
        return self._cache[key]

    def get_level(self, nuc: Nuclide) -> LevelData | None:
        """
        Decay data of a level, or None if the level is stable or has no data.

        Levels are matched as Geant4 does (G4VRadioactiveDecay::LoadDecayTable):
          * a level with a floating letter needs energy AND letter to match;
          * a level without a letter takes the first library level with the
            same energy, whatever its letter.  In that case the half-life is
            still the ENSDFSTATE one of the requested level, and transitions to
            the same nuclide are dropped (e.g. Tm164 ground state, 117 s, uses
            the decay modes of the library's '0 +X' level, 306 s, 80% IT);
          * an excited level (E > 0) absent from the library decays by IT to
            the ground state, with its ENSDFSTATE half-life;
          * a floating E = 0 level absent from the library (e.g. Ta178[0.000X],
            2.36 h) has no decay data: Geant4 kills it without emitting
            anything.  It is returned as None, i.e. terminal.
        """
        levels = self.get_levels(nuc.Z, nuc.A)
        for lv in levels:
            if abs(lv.excitation - nuc.E) < _E_TOL and lv.flb == nuc.flb:
                return lv
        if not nuc.flb:
            for lv in levels:
                if abs(lv.excitation - nuc.E) < _E_TOL:
                    hl = self.ensdf_halflife(nuc)
                    ds = [(d, br) for d, br in lv.daughters
                          if (d.Z, d.A) != (nuc.Z, nuc.A) or abs(d.E - nuc.E) >= _E_TOL]
                    norm = sum(br for _, br in ds)
                    self._warn(nuc, f"uses the decay data of library level "
                                    f"{Nuclide(nuc.Z, nuc.A, lv.excitation, lv.flb)}")
                    return lv._replace(
                        flb="",
                        half_life=hl if hl is not None else lv.half_life,
                        daughters=[(d, br / norm) for d, br in ds] if norm > 0 else [],
                    )
        hl = self.ensdf_halflife(nuc)
        if hl is None or hl <= 0:
            return None
        if nuc.E > 0:
            self._warn(nuc, "not in the decay library: IT to the ground state "
                            "(as in Geant4)", hl)
            gs = Nuclide(nuc.Z, nuc.A, 0.0, "")
            return LevelData(nuc.Z, nuc.A, nuc.E, nuc.flb, hl, [(gs, 1.0)], True)
        self._warn(nuc, "has no decay data: Geant4 kills it without emitting "
                        "radiation; treated as terminal", hl)
        return None

    def _warn(self, nuc: Nuclide, msg: str, hl: float | None = None) -> None:
        """Print a data warning once per level with 1 us <= T1/2 < 1e15 s (the
        upper limit skips double-beta emitters, stable for all purposes)."""
        if hl is None:
            hl = self.ensdf_halflife(nuc) or 0.0
        if self.verbose and 1e-6 <= hl < 1e15 and nuc not in self._warned:
            self._warned.add(nuc)
            print(f"  [WARNING] {nuc} (T1/2 = {hl:.4g} s) {msg}", file=sys.stderr)

    def __contains__(self, nuc: Nuclide) -> bool:
        return self.get_level(nuc) is not None


# ---------------------------------------------------------------------------
# Isotope name parser
# ---------------------------------------------------------------------------

def parse_isotope_name(name: str) -> Nuclide:
    """'Os172', 'Cs134[138.744]', 'Ta180[77.200X]' -> Nuclide (raises ValueError)."""
    return Nuclide(*parse_name(name))


# ---------------------------------------------------------------------------
# Core recursive chain builder
# ---------------------------------------------------------------------------

def _effective_daughters(
    nuc:          Nuclide,
    db:           DecayDatabase,
    min_halflife: float,
    ancestors:    frozenset[Nuclide],
) -> list[tuple[Nuclide, float]]:
    """
    The nuclides a decay into `nuc` effectively feeds: `nuc` itself if it is
    long-lived (T1/2 >= min_halflife) or stable; otherwise, recursively, the
    first long-lived or stable nuclides reached through its prompt decays,
    with their branching ratios.  Excited levels unknown to both data
    libraries collapse to the ground state.
    """
    gs = Nuclide(nuc.Z, nuc.A, 0.0, "")
    if nuc in ancestors:
        return [(nuc, 1.0)]
    level = db.get_level(nuc)
    if level is None:
        if nuc.E > 0 and db.ensdf_halflife(nuc) is None:
            # excited level unknown to both libraries: prompt de-excitation
            return _effective_daughters(gs, db, min_halflife, ancestors | {nuc})
        return [(nuc, 1.0)]                 # stable, or no decay data (terminal)
    if level.half_life < 0 or level.half_life >= min_halflife or not level.daughters:
        return [(nuc, 1.0)]
    out: dict[Nuclide, float] = {}
    for d, br in level.daughters:
        for e, br2 in _effective_daughters(d, db, min_halflife, ancestors | {nuc}):
            out[e] = out.get(e, 0.0) + br * br2
    return list(out.items())


def _recurse(
    nuc:             Nuclide,
    level:           LevelData,
    db:              DecayDatabase,
    cumulative_br:   float,
    current_path:    Chain,
    ancestors:       frozenset[Nuclide],
    results:         list[tuple[Chain, float]],
    prune_threshold: float,
    min_halflife:    float,
) -> None:
    """
    Depth-first traversal of the decay tree below the unstable nuclide `nuc`.

    Every recorded step goes from `nuc` (half-life level.half_life) to an
    effective daughter (see _effective_daughters), so that the daughter of
    step i is always the parent of step i+1.

    current_path : decay steps accumulated so far (modified in place, restored)
    ancestors    : nuclides already on this path (loop guard)
    results      : output list of (chain_copy, cumulative_br)
    """
    path_set = ancestors | {nuc}
    for daughter, br in level.daughters:
        for eff, br2 in _effective_daughters(daughter, db, min_halflife, path_set):
            step_br = br * br2
            cum = cumulative_br * step_br
            if cum < prune_threshold:
                continue
            current_path.append(DecayStep(level.half_life, eff, step_br))
            eff_level = None if eff in path_set else db.get_level(eff)
            if eff_level is None or eff_level.half_life < 0 or not eff_level.daughters:
                results.append((list(current_path), cum))      # terminal (or decay loop)
            else:
                _recurse(eff, eff_level, db, cum, current_path, path_set,
                         results, prune_threshold, min_halflife)
            current_path.pop()


def build_chain(
    isotope_name:    str,
    db:              DecayDatabase,
    prune_threshold: float = 0.001,
    min_halflife:    float = 1e-6,
) -> list[tuple[Chain, float]]:
    """
    Build all decay chains from a root isotope.

    Returns a list of (chain, cumulative_branching_ratio), sorted by BR
    descending and normalised to a total of 1.  Each chain is a list of
    DecayStep from the root to a terminal nuclide; the root itself is never
    skipped, whatever its half-life.  An empty list means the isotope is
    stable or absent from the database.
    """
    root = parse_isotope_name(isotope_name)
    level = db.get_level(root)
    if level is None or level.half_life < 0 or not level.daughters:
        return []

    results: list[tuple[Chain, float]] = []
    _recurse(root, level, db, 1.0, [], frozenset(), results, prune_threshold, min_halflife)

    # Merge chains with identical step sequences, summing their BRs.
    merged: dict[tuple, tuple[Chain, float]] = {}
    for chain, br in results:
        key = tuple((s.parent_halflife, s.daughter) for s in chain)
        if key in merged:
            merged[key] = (merged[key][0], merged[key][1] + br)
        else:
            merged[key] = (chain, br)
    results = [(c, br) for c, br in merged.values() if br >= prune_threshold]

    total = sum(br for _, br in results)
    if total > 0:
        results = [(c, br / total) for c, br in results]
    results.sort(key=lambda x: x[1], reverse=True)
    return results


# ---------------------------------------------------------------------------
# Output serialisation
# ---------------------------------------------------------------------------

def chain_filename(isotope_name: str) -> str:
    """Output filename for an isotope: 'Cs134[138.744]' -> 'Cs134[138.744].dat'."""
    return isotope_name + ".dat"


def params_line(db: DecayDatabase, prune_threshold: float, min_halflife: float) -> str:
    return (f"# params: format={CHAIN_FORMAT} prune={prune_threshold:g} "
            f"min_halflife={min_halflife:g} {db.tag}")


def read_params_line(path: str | Path) -> str | None:
    """The '# params:' line of an existing chain file, or None."""
    with open(path) as f:
        for line in f:
            if line.startswith("# params:"):
                return line.rstrip("\n")
            if not line.startswith("#"):
                break
    return None


def write_chain(
    isotope_name: str,
    chains:       list[tuple[Chain, float]],
    outdir:       str | Path = "DecayChains",
    params:       str | None = None,
) -> Path:
    """Write the chain list to <outdir>/<isotope>.dat (raises on empty chains)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not chains:
        raise ValueError(f"{isotope_name} has no decay chains (stable or not in database).")

    outpath = outdir / chain_filename(isotope_name)
    with outpath.open("w") as f:
        f.write(f"# {isotope_name} Decay Chain\n")
        if params:
            f.write(params + "\n")
        f.write("# HalfLife[s] (Z_daughter, A_daughter, E_daughter[, flb])\n")
        for i, (chain, br) in enumerate(chains):
            f.write(f"# Chain n. {i} Branching ratio: {br:.6f}\n")
            for step in chain:
                d = step.daughter
                flb = f", {d.flb}" if d.flb else ""
                f.write(f"{step.parent_halflife:.6e} ({d.Z}, {d.A}, {d.E}{flb})\n")
    return outpath


_STEP_RE = re.compile(
    r'([\deE+\-.]+)\s+\((\d+),\s*(\d+),\s*([\d.eE+\-]+)(?:,\s*([A-Z]))?\)')


def read_chain(path: str | Path) -> dict[int, tuple[float, list[tuple[float, Nuclide]]]]:
    """
    Read a chain file written by write_chain().

    Returns {chain_index: (branching_ratio, [(parent_halflife, daughter), ...])}.
    """
    chains: dict[int, tuple[float, list]] = {}
    current_idx: int | None = None
    current_steps: list = []
    current_br = 0.0

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or (line.startswith('#') and 'Chain n.' not in line):
            continue
        if 'Chain n.' in line:
            if current_idx is not None:
                chains[current_idx] = (current_br, current_steps)
            m = re.search(r'Chain n\.\s*(\d+)\s+Branching ratio:\s*([\d.eE+\-]+)', line)
            if m:
                current_idx, current_br, current_steps = int(m.group(1)), float(m.group(2)), []
        elif current_idx is not None:
            m = _STEP_RE.match(line)
            if m:
                nuc = Nuclide(int(m.group(2)), int(m.group(3)), float(m.group(4)),
                              m.group(5) or "")
                current_steps.append((float(m.group(1)), nuc))

    if current_idx is not None:
        chains[current_idx] = (current_br, current_steps)
    return chains


# ---------------------------------------------------------------------------
# Convenience: process a list of isotopes
# ---------------------------------------------------------------------------

def build_all(
    isotope_names:   list[str],
    db:              DecayDatabase,
    outdir:          str | Path = "DecayChains",
    prune_threshold: float = 0.001,
    min_halflife:    float = 1e-6,
    overwrite:       bool  = False,
    verbose:         bool  = True,
) -> dict[str, list[tuple[Chain, float]]]:
    """
    Build and save decay chains for a list of isotopes.

    An existing file is kept only if it was built with the same parameters
    and data libraries (its '# params:' line), unless overwrite=True.
    """
    outdir = Path(outdir)
    params = params_line(db, prune_threshold, min_halflife)
    results = {}
    for name in isotope_names:
        outpath = outdir / chain_filename(name)
        if outpath.exists() and not overwrite and read_params_line(outpath) == params:
            if verbose:
                print(f"  {name}: up to date, skipped.")
            continue
        if verbose:
            print(f"  {name}: building... ", end="", flush=True)
        try:
            chains = build_chain(name, db, prune_threshold=prune_threshold,
                                 min_halflife=min_halflife)
            if chains:
                write_chain(name, chains, outdir, params)
                results[name] = chains
                if verbose:
                    print(f"{len(chains)} chain(s) written to {outpath}")
            else:
                outpath.unlink(missing_ok=True)     # stale file from an older build
                if verbose:
                    print("stable or no decay data — no chain.")
        except Exception as exc:
            print(f"\n  [ERROR] {name}: {exc}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build radioactive decay chains from the Geant4 data libraries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("isotopes", nargs="+",
                        help="Isotope name(s), e.g. Os172 Pm139 Cs134[138.744]")
    parser.add_argument("--db", default=None,
                        help="RadioactiveDecay directory (default: $G4RADIOACTIVEDATA, "
                             "then $GEANT4_DATA_DIR/RadioactiveDecay*, then ./)")
    parser.add_argument("--ensdf", default=None,
                        help="G4ENSDFSTATE directory or ENSDFSTATE.dat (default: "
                             "$G4ENSDFSTATEDATA, then $GEANT4_DATA_DIR/G4ENSDFSTATE*)")
    parser.add_argument("--outdir", default="DecayChains", help="Output directory.")
    parser.add_argument("--prune", type=float, default=0.001, metavar="THRESHOLD",
                        help="Discard branches with cumulative branching ratio below this.")
    parser.add_argument("--min-halflife", type=float, default=1e-6, metavar="SECONDS",
                        help="Skip decay steps with half-life below this [s].")
    parser.add_argument("--overwrite", action="store_true",
                        help="Rebuild files even if up to date.")
    args = parser.parse_args(argv)

    db = DecayDatabase(args.db, args.ensdf)
    print(f"Decay library: {db._db_path}  ({db.size} nuclide files)")
    print(f"Half-lives:    {db._ensdf_path}\n")
    build_all(args.isotopes, db, outdir=args.outdir, prune_threshold=args.prune,
              min_halflife=args.min_halflife, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
