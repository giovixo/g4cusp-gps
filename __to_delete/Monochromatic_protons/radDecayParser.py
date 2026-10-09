import numpy as np
import re


"""
radDecayParser
Library for the parsing of radioactive data from Geant4 simulations

Version 1.0 - 2023-12-04

(c) Riccardo Campana - INAF/OAS
    riccardo.campana@inaf.it
"""


def split_element_string(inpstr):
    """
    Helper function to parse a Geant4 element string.
    
    Input:
        inpstr (string): input element string 
            (example: Gd145[749.100])
    
    Output:
        element, Z, A, excitationEnergy (tuple): output values 
            (example: Gd 64 143 749.100 for Gd-145 (Z=64) with 749.1 keV excitation)
    """
    periodic_table = {
            'H': 1,
            'He': 2,
            'Li': 3,
            'Be': 4,
            'B': 5,
            'C': 6,
            'N': 7,
            'O': 8,
            'F': 9,
            'Ne': 10,
            'Na': 11,
            'Mg': 12,
            'Al': 13,
            'Si': 14,
            'P': 15,
            'S': 16,
            'Cl': 17,
            'Ar': 18,
            'K': 19,
            'Ca': 20,
            'Sc': 21,
            'Ti': 22,
            'V': 23,
            'Cr': 24,
            'Mn': 25,
            'Fe': 26,
            'Co': 27,
            'Ni': 28,
            'Cu': 29,
            'Zn': 30,
            'Ga': 31,
            'Ge': 32,
            'As': 33,
            'Se': 34,
            'Br': 35,
            'Kr': 36,
            'Rb': 37,
            'Sr': 38,
            'Y': 39,
            'Zr': 40,
            'Nb': 41,
            'Mo': 42,
            'Tc': 43,
            'Ru': 44,
            'Rh': 45,
            'Pd': 46,
            'Ag': 47,
            'Cd': 48,
            'In': 49,
            'Sn': 50,
            'Sb': 51,
            'Te': 52,
            'I': 53,
            'Xe': 54,
            'Cs': 55,
            'Ba': 56,
            'La': 57,
            'Ce': 58,
            'Pr': 59,
            'Nd': 60,
            'Pm': 61,
            'Sm': 62,
            'Eu': 63,
            'Gd': 64,
            'Tb': 65,
            'Dy': 66,
            'Ho': 67,
            'Er': 68,
            'Tm': 69,
            'Yb': 70,
            'Lu': 71,
            'Hf': 72,
            'Ta': 73,
            'W': 74,
            'Re': 75,
            'Os': 76,
            'Ir': 77,
            'Pt': 78,
            'Au': 79,
            'Hg': 80,
            'Tl': 81,
            'Pb': 82,
            'Bi': 83,
            'Po': 84,
            'At': 85,
            'Rn': 86,
            'Fr': 87,
            'Ra': 88,
            'Ac': 89,
            'Th': 90,
            'Pa': 91,
            'U': 92,
            'Np': 93,
            'Pu': 94,
            'Am': 95,
            'Cm': 96,
            'Bk': 97,
            'Cf': 98,
            'Es': 99,
            'Fm': 100}

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
            excitationEnergy = 0
    else:
        print("ERROR!")
    for char in elemstr:
        if char.isalpha():
            alphabetic += char
        elif char.isnumeric():
            numeric_A += char
    return alphabetic, periodic_table[alphabetic], int(numeric_A), excitationEnergy
        

def formatElementString(Z, A, E):
    """
    Helper function to get a Geant4 element string from atomic and mass number
    
    Input:
        Z (int) = atomic number (example: 64)
        A (int) = mass number (example: 143)
        E (float) = excitation level in keV (example: 749.1)
    
    Output:
        output (string): output values (example: Gd143[749.1])
    """
    periodic_table = {
            'H': 1,
            'He': 2,
            'Li': 3,
            'Be': 4,
            'B': 5,
            'C': 6,
            'N': 7,
            'O': 8,
            'F': 9,
            'Ne': 10,
            'Na': 11,
            'Mg': 12,
            'Al': 13,
            'Si': 14,
            'P': 15,
            'S': 16,
            'Cl': 17,
            'Ar': 18,
            'K': 19,
            'Ca': 20,
            'Sc': 21,
            'Ti': 22,
            'V': 23,
            'Cr': 24,
            'Mn': 25,
            'Fe': 26,
            'Co': 27,
            'Ni': 28,
            'Cu': 29,
            'Zn': 30,
            'Ga': 31,
            'Ge': 32,
            'As': 33,
            'Se': 34,
            'Br': 35,
            'Kr': 36,
            'Rb': 37,
            'Sr': 38,
            'Y': 39,
            'Zr': 40,
            'Nb': 41,
            'Mo': 42,
            'Tc': 43,
            'Ru': 44,
            'Rh': 45,
            'Pd': 46,
            'Ag': 47,
            'Cd': 48,
            'In': 49,
            'Sn': 50,
            'Sb': 51,
            'Te': 52,
            'I': 53,
            'Xe': 54,
            'Cs': 55,
            'Ba': 56,
            'La': 57,
            'Ce': 58,
            'Pr': 59,
            'Nd': 60,
            'Pm': 61,
            'Sm': 62,
            'Eu': 63,
            'Gd': 64,
            'Tb': 65,
            'Dy': 66,
            'Ho': 67,
            'Er': 68,
            'Tm': 69,
            'Yb': 70,
            'Lu': 71,
            'Hf': 72,
            'Ta': 73,
            'W': 74,
            'Re': 75,
            'Os': 76,
            'Ir': 77,
            'Pt': 78,
            'Au': 79,
            'Hg': 80,
            'Tl': 81,
            'Pb': 82,
            'Bi': 83,
            'Po': 84,
            'At': 85,
            'Rn': 86,
            'Fr': 87,
            'Ra': 88,
            'Ac': 89,
            'Th': 90,
            'Pa': 91,
            'U': 92,
            'Np': 93,
            'Pu': 94,
            'Am': 95,
            'Cm': 96,
            'Bk': 97,
            'Cf': 98,
            'Es': 99,
            'Fm': 100}

    inv_periodic_table = {v: k for k, v in periodic_table.items()}
    if E > 0.:
        output = inv_periodic_table[Z] + "{:d}[{:.1f}]".format(A,E)
    else:
        output = inv_periodic_table[Z] + "{:d}".format(A)
    return output


def reformatElementName(s):
    """
    Helper function to reformat an element name
    
    Input:
        s (string) = input string (example: TODO)
    Output: 
        output (string) = output string (example: TODO)
    """
    A = [x for x in re.split('(\d+)',s) if x.isnumeric()][0]
    E = [x for x in re.split('(\d+)',s) if x.isalpha()][0]
    return E.title() + A
    

class ExcitationLevel(object):
    """
    Class which defines a radioactive isotope with its decay properties and excitation level
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
        
        periodic_table = {
            'H': 1,
            'He': 2,
            'Li': 3,
            'Be': 4,
            'B': 5,
            'C': 6,
            'N': 7,
            'O': 8,
            'F': 9,
            'Ne': 10,
            'Na': 11,
            'Mg': 12,
            'Al': 13,
            'Si': 14,
            'P': 15,
            'S': 16,
            'Cl': 17,
            'Ar': 18,
            'K': 19,
            'Ca': 20,
            'Sc': 21,
            'Ti': 22,
            'V': 23,
            'Cr': 24,
            'Mn': 25,
            'Fe': 26,
            'Co': 27,
            'Ni': 28,
            'Cu': 29,
            'Zn': 30,
            'Ga': 31,
            'Ge': 32,
            'As': 33,
            'Se': 34,
            'Br': 35,
            'Kr': 36,
            'Rb': 37,
            'Sr': 38,
            'Y': 39,
            'Zr': 40,
            'Nb': 41,
            'Mo': 42,
            'Tc': 43,
            'Ru': 44,
            'Rh': 45,
            'Pd': 46,
            'Ag': 47,
            'Cd': 48,
            'In': 49,
            'Sn': 50,
            'Sb': 51,
            'Te': 52,
            'I': 53,
            'Xe': 54,
            'Cs': 55,
            'Ba': 56,
            'La': 57,
            'Ce': 58,
            'Pr': 59,
            'Nd': 60,
            'Pm': 61,
            'Sm': 62,
            'Eu': 63,
            'Gd': 64,
            'Tb': 65,
            'Dy': 66,
            'Ho': 67,
            'Er': 68,
            'Tm': 69,
            'Yb': 70,
            'Lu': 71,
            'Hf': 72,
            'Ta': 73,
            'W': 74,
            'Re': 75,
            'Os': 76,
            'Ir': 77,
            'Pt': 78,
            'Au': 79,
            'Hg': 80,
            'Tl': 81,
            'Pb': 82,
            'Bi': 83,
            'Po': 84,
            'At': 85,
            'Rn': 86,
            'Fr': 87,
            'Ra': 88,
            'Ac': 89,
            'Th': 90,
            'Pa': 91,
            'U': 92,
            'Np': 93,
            'Pu': 94,
            'Am': 95,
            'Cm': 96,
            'Bk': 97,
            'Cf': 98,
            'Es': 99,
            'Fm': 100
        }
        self.isoA =  int(''.join([i for i in self.name if i.isdigit()]))
        element_name = ''.join([i for i in self.name if not i.isdigit()])
        self.isoZ = int(periodic_table[element_name])
        
    def __str__(self):
        str_out = "Element: {:s} ({:d}, {:d}) \t Excitation level: {:.2f} keV \t Half life: {:.3e} s\n".format(self.name, self.isoZ, self.isoA, self.excitationEnergy, self.halfLife)
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
        newDecayModes = []
        for d in self.decayModes:
            if len(d) == 3:
                mode, daughterEx, intensity = d
                newDecayModes.append([mode, float(daughterEx), float(intensity)])
            elif len(d) == 4:
                mode, daughterEx, flag, intensity = d
                newDecayModes.append([mode, float(daughterEx), float(intensity)])
            #elif len(d) == 5:
                #mode, daughterEx, flag, intensity, q = d
            #elif len(d) == 6:
                #mode, daughterEx, flag, intensity, q, note = d
            elif len(d)>7:
                print("ANOMALOUS LENGTH", len(d))
        s = 0        
        for nd in newDecayModes:
            s += nd[2]
        for nd in newDecayModes:
            # Renormalise intensity to avoid rounding errors in branching ratios
            nd[2] = nd[2]/s
        self.decayModes = newDecayModes
        
    def identifyDaughters(self):
        daughters = []
        
        for decaymode in self.decayModes:
            if decaymode[0] == "Alpha":
                # Alpha decay
                daughterZ = self.isoZ-2
                daughterA = self.isoA-4
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "BetaPlus" or decaymode[0][-2:] == "EC":
                # Beta+ decay or electronic capture
                daughterZ = self.isoZ-1
                daughterA = self.isoA
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "BetaMinus":
                # Beta- decay 
                daughterZ = self.isoZ+1
                daughterA = self.isoA
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "IT":
                # Isomeric transition
                daughterZ = self.isoZ
                daughterA = self.isoA
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "SpFission":
                # Spontaneous fission: skipped
                print("Decay mode SpFission --> SKIPPED")
            elif decaymode[0] == "Neutron":
                daughterZ = self.isoZ
                daughterA = self.isoA-1
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "Proton":
                daughterZ = self.isoZ-1
                daughterA = self.isoA-1
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "Beta2Minus":
                daughterZ = self.isoZ+2
                daughterA = self.isoA
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "Beta2Plus":
                daughterZ = self.isoZ-2
                daughterA = self.isoA
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            elif decaymode[0] == "Triton":
                daughterZ = self.isoZ
                daughterA = self.isoA-3
                daughters.append([(daughterZ, daughterA, decaymode[1]), decaymode[2]])
            else:
                # TODO: Put here other decay modes if any
                print("Decay mode", decaymode[0], "--> SKIPPED")
        
        # Consolidate states that will decay in the same daughter isotope
        # identified by unique Z,A,E

        set_daughtersZAE = set([d[0] for d in daughters]) 
        newDaughters = []
        for k in set_daughtersZAE:
            totChannelIntensity = 0
            for d in daughters:
                if d[0] == k:
                    totChannelIntensity += d[1]
            newDaughters.append([k, self.halfLife, totChannelIntensity])
        
        if len(newDaughters) == 1:
           newDaughters[0][2] = 1.
        self.daughters = newDaughters
        
    def sumIntensities(self):
        s = 0
        for d in self.decayModes:
            if len(d) == 3:
                mode, daughterEx, intensity = d
            elif len(d) == 5:
                mode, daughterEx, flag, intensity, q = d
            elif len(d) == 6:
                mode, daughterEx, flag, intensity, q, note = d
            else:
                print("ANOMALOUS LENGTH", len(d))
            s += float(intensity)
        return(s)



def getDecay(Z, A, excitationEnergy=0., verbose=False):
    """
    Parser function for the Geant4 RadioactiveDecay database
    
    Input:
        Z (int) = Atomic number
        A (int) = Mass number
        excitationEnergy (float) = Excited level energy in keV
    
    Output:
        output (ExcitationLevel object) = ExcitationLevel() class object corresponding to the decay
    
    TODO: treat floating levels, now they are skipped.
    """
    
    
    filein = "./RadioactiveDecay6.1.2/z{:d}.a{:d}".format(Z, A)

    f = open(filein)
    lines = f.readlines()
    f.close()
    
    levels = []
    exlevel = None
    decays = None

    for i,line in enumerate(lines):
        if i == 0 and line[0] == '#':
            # First comment line
            # print(line)
            d = line.split()
            # print(d)
            elementName = reformatElementName(d[1])
        elif line[0] == '#':
            pass
        elif line[0] == 'P':
            # Found first line of excitation level entry
            if exlevel is not None:
                exlevel.decayModes = decays
                levels.append(exlevel)
            exlevel = None
            del decays
            p = line.split()
            if len(lines) < 4:
                print("Wtf is happening there??!")
                print(Z, A)
            excitationLevel = float(p[1])
            flag = p[2]
            halfLife = float(p[3])
            if halfLife < 0:
                print("HALF LIFE NEGATIVE!!! Decay chain stops here.")
                skip = True
            else:
                if flag != '-':
                    exlevel = ExcitationLevel(elementName, excitationLevel, halfLife, floating=flag)
                else:
                    exlevel = ExcitationLevel(elementName, excitationLevel, halfLife)
                skip = False
            #print(exlevel)
            decays = []
        else:
            #print(line.split())
            if not skip:
                decays.append(line.split())
    if exlevel is not None:
        exlevel.decayModes = decays
        
    levels.append(exlevel)

    output = None

    if levels:
        for exl in levels:
            exl.reformatDecayModes()
            exl.identifyDaughters()
        
    if verbose:
        for exl in levels:
            print(exl)

    for exl in levels:
        if (exl.excitationEnergy-excitationEnergy < 1e-3) and (not exl.isFloating):
            if verbose:
                print("Found an entry with excitation energy:", exl.excitationEnergy)
                print("which is:")
                print(exl)
            output = exl
    return output
    
# print()
# exl = getDecay(65,149, excitationEnergy=0)
# print(exl)

def get_decay_constant(isotopeName):
    iso, Z, A, E = split_element_string(isotopeName)
    decay = getDecay(Z, A, excitationEnergy=E)
    try:
        halflife = decay.halfLife
        l = 0.693/halflife
        return l
    except AttributeError:
        print("Tb143/Lu162 bug?")
        print(isotopeName)
        return 0.
    except IndexError:
        print("Tb143/Lu162 bug?")
        print(isotopeName)
        return 0.
        
def get_daughters(isotopeName):
    try:
        f = open("DecayChains/"+isotopeName+".dat", "r")
        lines = f.readlines()
        f.close()
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
            isoE = float(line.split("(")[1].split(",")[2].replace(")",""))
            d[br].append([half_life, (isoZ, isoA, isoE)])
    return d
    
    # if len(lines) > 3:
    #     daughters = [x.strip() for x in lines[3:]]
    #     print("Decay chain:")
    #     for daughter in daughters:
    #         print(daughter)
    #         half_life = float(daughter.split("(")[0])
    #         isoZ = int(daughter.split("(")[1].split(",")[0])
    #         isoA = int(daughter.split("(")[1].split(",")[1])
    #         isoE = float(daughter.split("(")[1].split(",")[2].replace(")",""))
    #         print(half_life, isoZ, isoA, isoE)