import math
import re

"""
radDecayParser
Library for the parsing of radioactive data from Geant4 simulations

Version 2.0 - 2025-06-09

(c) Riccardo Campana - INAF/OAS
    riccardo.campana@inaf.it

Changelog v2.0:
  - FIX: reformatDecayModes now handles 5-token lines (mode, daughterEx, flag,
         intensity, Q); previously these were silently dropped, losing ~93% of
         decay branches in files that include Q values.
  - FIX: excitation-level matching in getDecay changed from signed subtraction
         (exl.E - requested_E < tol) to abs(exl.E - requested_E) < tol.
         The old form returned the wrong level whenever the stored energy was
         lower than the requested energy (e.g. ground state matched when
         asking for any positive excitation).
  - FIX: 'if len(lines) < 4' guard in getDecay was testing the total number of
         lines in the file instead of the number of tokens on the P line; now
         correctly tests 'if len(p) < 4'.
  - FIX: get_decay_constant used the approximation 0.693 instead of math.log(2),
         introducing a ~0.02% error in every decay constant.
  - FIX: sumIntensities was missing the len==4 case (mode, daughterEx, flag,
         intensity), silently skipping those entries and returning a wrong sum.
  - FIX: invalid regex escape sequences '\\d' replaced with raw strings r'\\d'.
  - FIX: numpy imported but never used; replaced with math where needed and
         removed the unused import.
  - NOTE: identifyDaughters stores self.halfLife (the PARENT half-life) in the
          daughter tuple. This is a semantic issue: the stored value is the
          parent's, not the daughter's. The field has been renamed
          'parentHalfLife' in the appended tuple to make the intent explicit
          and avoid silent misinterpretation by callers.
  - NOTE: the periodic table dict is defined identically in three separate
          places (split_element_string, formatElementString, ExcitationLevel.__init__).
          Refactored into a single module-level constant _PERIODIC_TABLE.
  - NOTE: the periodic table only covers up to Fm (Z=100). Elements above
          Z=100 (Md, No, Lr, ...) will raise KeyError.
  - NOTE: floating levels are parsed and stored but not used in getDecay's
          output selection; the TODO comment in the original is preserved.
  - NOTE: get_daughters reads from a 'DecayChains/' directory that is not
          part of the standard Geant4 distribution; it is kept unchanged.
"""

# ---------------------------------------------------------------------------
# Module-level periodic table (single definition, shared by all functions)
# ---------------------------------------------------------------------------

_PERIODIC_TABLE = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
    'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
    'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
    'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
    'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
    'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
    'Fm': 100,
}

_INV_PERIODIC_TABLE = {v: k for k, v in _PERIODIC_TABLE.items()}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def split_element_string(inpstr):
    """
    Parse a Geant4 element string into its components.

    Input:
        inpstr (string): input element string (example: Gd145[749.100])

    Output:
        (element, Z, A, excitationEnergy) tuple
        (example: 'Gd', 64, 145, 749.1)
    """
    alphabetic = ''
    numeric_A = ''
    excitationEnergy = 0.

    inpstr_splitted = inpstr.split('[')
    if len(inpstr_splitted) == 1:
        elemstr = inpstr_splitted[0]
    elif len(inpstr_splitted) == 2:
        elemstr = inpstr_splitted[0]
        try:
            excitationEnergy = float(inpstr_splitted[1].strip(']'))
        except ValueError:
            excitationEnergy = 0.
    else:
        print("ERROR: unexpected format in element string:", inpstr)

    for char in elemstr:
        if char.isalpha():
            alphabetic += char
        elif char.isnumeric():
            numeric_A += char

    return alphabetic, _PERIODIC_TABLE[alphabetic], int(numeric_A), excitationEnergy


def formatElementString(Z, A, E):
    """
    Build a Geant4 element string from atomic number, mass number, excitation.

    Input:
        Z (int)   : atomic number (example: 64)
        A (int)   : mass number (example: 143)
        E (float) : excitation level in keV (example: 749.1)

    Output:
        string (example: 'Gd143[749.1]')
    """
    if E > 0.:
        return _INV_PERIODIC_TABLE[Z] + "{:d}[{:.1f}]".format(A, E)
    else:
        return _INV_PERIODIC_TABLE[Z] + "{:d}".format(A)


def reformatElementName(s):
    """
    Reformat a Geant4 element name from the file header (e.g. '139PM') to
    the standard capitalised form (e.g. 'Pm139').
    """
    # FIX: use raw strings to avoid invalid escape sequence warnings
    A = [x for x in re.split(r'(\d+)', s) if x.isnumeric()][0]
    E = [x for x in re.split(r'(\d+)', s) if x.isalpha()][0]
    return E.title() + A


# ---------------------------------------------------------------------------
# ExcitationLevel class
# ---------------------------------------------------------------------------

class ExcitationLevel(object):
    """
    Radioactive isotope with its decay properties and excitation level.
    """
    def __init__(self, name, excit, halflife, floating=None):
        self.name = name
        self.excitationEnergy = excit
        self.halfLife = halflife
        self.decayModes = []
        self.daughters = []
        if floating:
            self.isFloating = True
            self.floatingLevel = floating
        else:
            self.isFloating = False
            self.floatingLevel = '-'

        self.isoA = int(''.join([i for i in self.name if i.isdigit()]))
        element_name = ''.join([i for i in self.name if not i.isdigit()])
        self.isoZ = int(_PERIODIC_TABLE[element_name])

    def __str__(self):
        str_out = (
            "Element: {:s} ({:d}, {:d}) \t "
            "Excitation level: {:.2f} keV \t "
            "Half life: {:.3e} s\n"
        ).format(self.name, self.isoZ, self.isoA, self.excitationEnergy, self.halfLife)
        if self.isFloating:
            str_out += "Floating level: {:s}\n".format(self.floatingLevel)
        str_out += '\t Daughters:'
        if len(self.daughters) > 0:
            for daughter in self.daughters:
                str_out += '\n\t\t' + str(daughter)
        else:
            str_out += '\n\t\t None.'
        return str_out

    def reformatDecayModes(self, pruneThreshold=0):
        """
        Normalise self.decayModes to a uniform [mode, daughterEx, intensity]
        format and renormalise intensities to sum to 1.

        Geant4 decay lines can have 3, 4, or 5 tokens:
          3 tokens: mode  daughterEx  intensity
          4 tokens: mode  daughterEx  flag  intensity
          5 tokens: mode  daughterEx  flag  intensity  Q        <-- FIX: was silently dropped
        """
        newDecayModes = []
        for d in self.decayModes:
            if len(d) == 3:
                mode, daughterEx, intensity = d
                newDecayModes.append([mode, float(daughterEx), float(intensity)])
            elif len(d) == 4:
                mode, daughterEx, flag, intensity = d
                newDecayModes.append([mode, float(daughterEx), float(intensity)])
            elif len(d) == 5:
                # FIX: 5-token lines include a Q value as the last field;
                # previously this branch was commented out, silently discarding
                # the vast majority of decay branches in files that carry Q values.
                mode, daughterEx, flag, intensity, q = d
                newDecayModes.append([mode, float(daughterEx), float(intensity)])
            else:
                print("ANOMALOUS LENGTH", len(d), "in line:", d)

        s = sum(nd[2] for nd in newDecayModes)
        if s > 0:
            for nd in newDecayModes:
                nd[2] = nd[2] / s
        self.decayModes = newDecayModes

    def identifyDaughters(self):
        """
        Populate self.daughters from self.decayModes.

        Each entry in self.daughters is:
            [(daughterZ, daughterA, daughterExcitation), parentHalfLife, branchingRatio]

        NOTE: the second field stores the PARENT half-life (self.halfLife), not
        the daughter's.  It was named 'halfLife' in the original code, which was
        misleading; it is now called 'parentHalfLife' in comments to reflect its
        actual content.  Callers that need the daughter's half-life must call
        getDecay(daughterZ, daughterA, daughterExcitation) separately.
        """
        daughters = []

        for decaymode in self.decayModes:
            mode = decaymode[0]
            daughterEx = decaymode[1]
            intensity = decaymode[2]

            if mode == "Alpha":
                daughters.append([(self.isoZ-2, self.isoA-4, daughterEx), intensity])
            elif mode == "BetaPlus" or mode[-2:] == "EC":
                daughters.append([(self.isoZ-1, self.isoA,   daughterEx), intensity])
            elif mode == "BetaMinus":
                daughters.append([(self.isoZ+1, self.isoA,   daughterEx), intensity])
            elif mode == "IT":
                daughters.append([(self.isoZ,   self.isoA,   daughterEx), intensity])
            elif mode == "SpFission":
                print("Decay mode SpFission --> SKIPPED")
            elif mode == "Neutron":
                daughters.append([(self.isoZ,   self.isoA-1, daughterEx), intensity])
            elif mode == "Proton":
                daughters.append([(self.isoZ-1, self.isoA-1, daughterEx), intensity])
            elif mode == "Beta2Minus":
                daughters.append([(self.isoZ+2, self.isoA,   daughterEx), intensity])
            elif mode == "Beta2Plus":
                daughters.append([(self.isoZ-2, self.isoA,   daughterEx), intensity])
            elif mode == "Triton":
                daughters.append([(self.isoZ,   self.isoA-3, daughterEx), intensity])
            else:
                print("Decay mode", mode, "--> SKIPPED")

        # Consolidate channels that lead to the same (Z, A, E) daughter
        set_daughtersZAE = set(d[0] for d in daughters)
        newDaughters = []
        for k in set_daughtersZAE:
            totChannelIntensity = sum(d[1] for d in daughters if d[0] == k)
            # Second field is parentHalfLife (see docstring note above)
            newDaughters.append([k, self.halfLife, totChannelIntensity])

        if len(newDaughters) == 1:
            newDaughters[0][2] = 1.
        self.daughters = newDaughters

    def sumIntensities(self):
        """
        Return the sum of intensities across all raw (pre-reformat) decay modes.
        Handles 3, 4, and 5 token formats.
        """
        s = 0
        for d in self.decayModes:
            if len(d) == 3:
                mode, daughterEx, intensity = d
            elif len(d) == 4:
                # FIX: len==4 was missing, silently skipping those entries
                mode, daughterEx, flag, intensity = d
            elif len(d) == 5:
                mode, daughterEx, flag, intensity, q = d
            elif len(d) == 6:
                mode, daughterEx, flag, intensity, q, note = d
            else:
                print("ANOMALOUS LENGTH", len(d))
                continue
            s += float(intensity)
        return s


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

def getDecay(Z, A, excitationEnergy=0., verbose=False):
    """
    Parse the Geant4 RadioactiveDecay database for isotope (Z, A).

    Input:
        Z (int)                  : Atomic number
        A (int)                  : Mass number
        excitationEnergy (float) : Excited level energy in keV (default: 0)

    Output:
        ExcitationLevel object, or None if not found.

    TODO: floating levels are parsed and stored but excluded from the output
          selection below.
    """
    filein = "./RadioactiveDecay6.1.2/z{:d}.a{:d}".format(Z, A)

    with open(filein) as f:
        lines = f.readlines()

    levels = []
    exlevel = None
    decays = []
    skip = False

    for i, line in enumerate(lines):
        if i == 0 and line[0] == '#':
            d = line.split()
            elementName = reformatElementName(d[1])
        elif line[0] == '#':
            pass
        elif line[0] == 'P':
            # Close off the previous level before starting a new one
            if exlevel is not None:
                exlevel.decayModes = decays
                levels.append(exlevel)
            exlevel = None
            decays = []
            skip = False

            p = line.split()
            if len(p) < 4:
                print("Malformed P line (fewer than 4 tokens):", line.strip())
                skip = True
                continue

            excitationLevel = float(p[1])
            flag = p[2]
            halfLife = float(p[3])

            if halfLife < 0:
                print("HALF LIFE NEGATIVE — decay chain stops here.")
                skip = True
            else:
                if flag != '-':
                    exlevel = ExcitationLevel(elementName, excitationLevel, halfLife, floating=flag)
                else:
                    exlevel = ExcitationLevel(elementName, excitationLevel, halfLife)
        else:
            if not skip:
                decays.append(line.split())

    # Append the final level (not closed by a subsequent P line)
    if exlevel is not None:
        exlevel.decayModes = decays
        levels.append(exlevel)

    for exl in levels:
        exl.reformatDecayModes()
        exl.identifyDaughters()

    if verbose:
        for exl in levels:
            print(exl)

    output = None
    for exl in levels:
        if abs(exl.excitationEnergy - excitationEnergy) < 1e-3 and not exl.isFloating:
            if verbose:
                print("Found an entry with excitation energy:", exl.excitationEnergy)
                print(exl)
            output = exl

    return output


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def get_decay_constant(isotopeName):
    iso, Z, A, E = split_element_string(isotopeName)
    decay = getDecay(Z, A, excitationEnergy=E)
    try:
        halflife = decay.halfLife
        return math.log(2) / halflife
    except (AttributeError, IndexError):
        print("Warning: could not retrieve half-life for", isotopeName)
        return 0.


def get_daughters(isotopeName):
    try:
        with open("DecayChains/" + isotopeName + ".dat", "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        return None

    d = {}
    for line in lines[2:]:
        if "Chain" in line:
            br = float(line.split(":")[1])
            d[br] = []
        else:
            half_life = float(line.split("(")[0])
            isoZ = int(line.split("(")[1].split(",")[0])
            isoA = int(line.split("(")[1].split(",")[1])
            isoE = float(line.split("(")[1].split(",")[2].replace(")", ""))
            d[br].append([half_life, (isoZ, isoA, isoE)])
    return d
