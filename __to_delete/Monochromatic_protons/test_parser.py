import os
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional


# ---------------------------------------------------------
# Nuclear bookkeeping
# ---------------------------------------------------------

ELEMENTS = {
    1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O",
    9:"F", 10:"Ne", 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P",
    16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca", 21:"Sc", 22:"Ti",
    23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu",
    30:"Zn", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr",
    37:"Rb", 38:"Sr", 39:"Y", 40:"Zr", 41:"Nb", 42:"Mo", 43:"Tc",
    44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn",
    51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La",
    58:"Ce", 59:"Pr", 60:"Nd", 61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd",
    65:"Tb", 66:"Dy", 67:"Ho", 68:"Er", 69:"Tm", 70:"Yb", 71:"Lu",
    72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt",
    79:"Au", 80:"Hg", 81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At",
    86:"Rn", 87:"Fr", 88:"Ra", 89:"Ac", 90:"Th", 91:"Pa", 92:"U",
    93:"Np", 94:"Pu"
}


# ---------------------------------------------------------
# Decay mode daughter calculation
# ---------------------------------------------------------

def daughter_nuclide(Z, A, mode):
    """
    Determine daughter nuclide from decay mode.
    """

    mode = mode.lower()

    if "betaminus" in mode:
        return Z + 1, A

    if "betaplus" in mode or "ec" in mode:
        return Z - 1, A

    if mode == "alpha":
        return Z - 2, A - 4

    if mode == "it":
        return Z, A

    if mode == "proton":
        return Z - 1, A - 1

    if mode == "neutron":
        return Z, A - 1

    return None


# ---------------------------------------------------------
# Data structures
# ---------------------------------------------------------

@dataclass
class Branch:
    mode: str
    intensity: float
    daughter_Z: int
    daughter_A: int


@dataclass
class Nuclide:
    Z: int
    A: int
    branches: List[Branch]


# ---------------------------------------------------------
# Parser
# ---------------------------------------------------------

class G4RadioactiveData:

    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.cache: Dict[tuple, Nuclide] = {}

    def filename(self, Z, A):
        return self.data_dir / f"z{Z}.a{A}"

    def load(self, Z, A) -> Optional[Nuclide]:

        key = (Z, A)

        if key in self.cache:
            return self.cache[key]

        fn = self.filename(Z, A)

        if not fn.exists():
            return None

        branches = []

        with open(fn) as f:
            for line in f:

                line = line.strip()

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                tokens = line.split()

                if len(tokens) == 0:
                    continue

                mode = tokens[0]

                # Look for decay records
                if mode in (
                    "BetaMinus",
                    "BetaPlus",
                    "Alpha",
                    "IT",
                    "EC",
                    "Proton",
                    "Neutron"
                ):

                    try:
                        intensity = float(tokens[-2])
                    except Exception:
                        intensity = 0.0

                    daughter = daughter_nuclide(Z, A, mode)

                    if daughter is None:
                        continue

                    dZ, dA = daughter

                    branches.append(
                        Branch(
                            mode=mode,
                            intensity=intensity,
                            daughter_Z=dZ,
                            daughter_A=dA
                        )
                    )

        nuclide = Nuclide(Z, A, branches)

        self.cache[key] = nuclide

        return nuclide


# ---------------------------------------------------------
# Recursive decay chain
# ---------------------------------------------------------

def build_chain(db, Z, A, visited=None):

    if visited is None:
        visited = set()

    key = (Z, A)

    if key in visited:
        return {
            "nuclide": key,
            "loop": True
        }

    visited.add(key)

    nuclide = db.load(Z, A)

    if nuclide is None:
        return {
            "nuclide": key,
            "stable_or_missing": True
        }

    node = {
        "nuclide": key,
        "branches": []
    }

    for b in nuclide.branches:

        child = build_chain(
            db,
            b.daughter_Z,
            b.daughter_A,
            visited.copy()
        )

        node["branches"].append(
            {
                "mode": b.mode,
                "intensity": b.intensity,
                "child": child
            }
        )

    return node


# ---------------------------------------------------------
# Pretty printer
# ---------------------------------------------------------

def nuclide_name(Z, A):
    return f"{A}{ELEMENTS.get(Z, f'Z{Z}')}"


def print_chain(node, indent=""):

    Z, A = node["nuclide"]

    print(indent + nuclide_name(Z, A))

    if node.get("stable_or_missing"):
        return

    for branch in node["branches"]:

        child = branch["child"]

        dZ, dA = child["nuclide"]

        print(
            indent +
            f"  └─[{branch['mode']}, I={branch['intensity']}]→ "
            f"{nuclide_name(dZ, dA)}"
        )

        print_chain(child, indent + "      ")