"""
nuclides.py
===========
Canonical nuclide names shared by all pipeline steps.

Geant4 names ions as  <Sym><A>  for ground states and  <Sym><A>[<E>]  for
excited states, with E in keV printed with 5 decimals and optionally followed
by a floating-level-base letter (G4Ions::FloatLevelBaseChar: X, Y, Z, U, V, W,
R, S, T, A, B, C, D, E), e.g. 'Al26[228.30500]' or 'Ta180[77.20000X]'.

The pipeline uses ONE canonical form everywhere (activation tables, decay
chain files, activity tables, spectrum file names):

    Ga70            ground state
    Al26[228.305]   excited state, E in keV with exactly 3 decimals
    Ta180[77.200X]  excited state on a floating level

Always build names with format_name() and read them with parse_name().
"""

from __future__ import annotations

import re
from typing import NamedTuple

PERIODIC_TABLE: dict[str, int] = {
    'H': 1,  'He': 2,  'Li': 3,  'Be': 4,  'B': 5,   'C': 6,   'N': 7,
    'O': 8,  'F': 9,   'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14,
    'P': 15, 'S': 16,  'Cl': 17, 'Ar': 18, 'K': 19,  'Ca': 20, 'Sc': 21,
    'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35,
    'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39,  'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
    'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53,  'Xe': 54, 'Cs': 55, 'Ba': 56,
    'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63,
    'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,  'Re': 75, 'Os': 76, 'Ir': 77,
    'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
    'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
    'U': 92,  'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
    'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104,
    'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110,
    'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
    'Ts': 117, 'Og': 118,
}
SYMBOLS: dict[int, str] = {z: s for s, z in PERIODIC_TABLE.items()}

FLOAT_LEVEL_CHARS = "XYZUVWRSTABCDE"

_NAME_RE = re.compile(
    r"^(?P<sym>[A-Z][a-z]?)(?P<A>\d+)"
    r"(?:\[(?P<E>\d+(?:\.\d*)?)(?P<flb>[" + FLOAT_LEVEL_CHARS + r"]?)\])?$"
)


class NuclideId(NamedTuple):
    Z: int
    A: int
    E: float        # excitation energy [keV]
    flb: str = ""   # floating level base letter ('' = none)


def format_name(Z: int, A: int, E: float = 0.0, flb: str = "") -> str:
    """Canonical name: 'Ga70', 'Al26[228.305]', 'Ta180[77.200X]'."""
    sym = SYMBOLS.get(Z, f"Z{Z}")
    if E > 0.0 or flb:
        return f"{sym}{A}[{E:.3f}{flb}]"
    return f"{sym}{A}"


def parse_name(name: str) -> NuclideId:
    """Parse a Geant4 or canonical nuclide name. Raises ValueError if malformed."""
    m = _NAME_RE.match(name.strip())
    if m is None:
        raise ValueError(f"Cannot parse nuclide name {name!r}")
    sym = m.group("sym")
    if sym not in PERIODIC_TABLE:
        raise ValueError(f"Unknown element symbol {sym!r} in {name!r}")
    E = float(m.group("E")) if m.group("E") else 0.0
    return NuclideId(PERIODIC_TABLE[sym], int(m.group("A")), E, m.group("flb") or "")


def canonical_name(name: str) -> str:
    """Geant4 name -> canonical name ('Al26[228.30500]' -> 'Al26[228.305]')."""
    return format_name(*parse_name(name))
