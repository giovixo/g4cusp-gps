# RadioactiveDecay6.1.3d: README_RDM does not describe the file format that the files and LoadDecayTable use

- **Component:** data library RadioactiveDecay6.1.3d (README_RDM) and processes/hadronic/models/radioactive_decay (`G4VRadioactiveDecay::LoadDecayTable`)
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2 and 6.1.3d, G4ENSDFSTATE3.0)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** minor. Misleading documentation; files edited by hand following it can be read wrongly by Geant4

README_RDM is the same in RadioactiveDecay6.1.2 and 6.1.3d. Two other problems found in the 6.1.3 files of 29 October 2026 are fixed in 6.1.3d: the 138-character comment line in `z86.a197`, which made `LoadDecayTable` fail with warning HAD_RDM_100 (6.1.3d: longest line 115 characters, no warning with rdecay01), and the half-lives in the comment headers that differed from ENSDFSTATE (6.1.3d: none of the 2209 headers differs by more than 1%).

## Summary

README_RDM describes a format that differs from the one used by the files and read by `LoadDecayTable`:

1. It gives the summary lines four columns; they have three.
2. It does not document the fifth column of the detail lines, the Q-value, which Geant4 uses.
3. It says that the detail percentages are relative to the mode total; in about half of the levels they are not.
4. It does not say that `LoadDecayTable` tells the line types apart by their length, nor that it cannot read lines of 120 characters or more.
5. It says that the half-life of the `P` lines is ignored only for some decay modes; it is ignored for all of them.
6. It does not list the decay-mode keywords.

The numbers below come from the attached script, for 6.1.3d (6.1.2 in parentheses).

## Details

1. **Summary lines.** README_RDM: "If there are four columns only, this is a header describing the total strength of each type of decay from this level. The first column is the decay mode, the second column is always zero, the third column indicates whether the level is floating or not, … and the fourth column is the total branching ratio". The files have three columns (mode, 0, fraction): 9861 (9832) three-column lines and no four-column ones, and no floating flag, either as a separate column or attached to another one. `LoadDecayTable` reads three fields (`tmpStream >> theDecayMode >> dummy >> decayModeTotal`, `G4VRadioactiveDecay.cc:627`).

2. **Q column of the detail lines.** README_RDM describes the detail lines as five columns but names only four fields (decay type, daughter level, floating-level indicator, branching ratio). The fifth column is the Q-value in keV. `LoadDecayTable` reads it (`tmpStream >> theDecayMode >> a >> daughterFloatFlag >> b >> c`, `G4VRadioactiveDecay.cc:679`) and passes it to the decay channel as the Q-value, so a wrong value changes the emitted energies or suppresses the decay (bug 2780). The optional sixth column, the β type, is documented.

3. **Normalisation of the detail percentages.** README_RDM says that the detail branching ratio is "a percentage … of the total branching ratio which came from the 4-column record for this decay", i.e. the detail lines of each mode add up to 100. In 1495 of the 3241 (1486 of 3228) levels with detail lines, the percentages of at least one mode do not add up to 100 (±1); they are relative to all decays of the level. This is harmless while the summary fractions are right, because `LoadDecayTable` rescales the detail lines of each mode to the summary fraction of that mode (but see report D).

4. **Line types and line length.** `LoadDecayTable` does not count columns. After removing trailing blanks, it takes a line shorter than 72 characters as a summary line, one shorter than 84 as a detail line without β type, and a longer one as a detail line with β type (`G4VRadioactiveDecay.cc:626`, `678`). It reads the file with `getline(inputChars, 120)` (`G4VRadioactiveDecay.cc:567`, `589`): a line of 120 characters or more, comments included, makes `getline` fail, the rest of the file is not read, and the loop ends only at its guard (warning HAD_RDM_100). This happened with a comment line of `z86.a197` in the 6.1.3 files of 29 October 2026. README_RDM mentions neither rule, so a file edited by hand by following it can easily break them.

5. **Half-life of the `P` lines.** README_RDM: "The number following this is the half life in seconds. This value is ignored for alpha and beta decay, and internal transitions. Instead, half life (mean life) is taken from the file ENSDFSTATE.dat". `LoadDecayTable` ignores it for every mode (`tmpStream >> recordType >> parentExcitation >> floatingFlag >> dummy`, with the comment "dummy takes the place of half-life", `G4VRadioactiveDecay.cc:605`).

6. **Decay-mode keywords.** README_RDM does not list the keywords of the first column. The files use 12: `IT`, `BetaMinus`, `BetaPlus`, `KshellEC`, `LshellEC`, `MshellEC`, `NshellEC`, `Alpha`, `Proton`, `Neutron`, `SpFission`, `Triton`. `G4RadioactiveDecayMode.cc` accepts six more (`BDProton`, `BDNeutron`, `Beta2Minus`, `Beta2Plus`, `Proton2`, `Neutron2`).

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_E_halflives_readme.py <RadioactiveDecay6.1.3d dir>   # default: 6.1.2 from $GEANT4_DATA_DIR
```

Section 2 of the output ("File format") gives the column counts, the line lengths and the levels whose detail percentages do not add up to 100. Section 1 checks the half-lives of the `P` lines and of the comment headers against ENSDFSTATE.

Attached:
- `check_E_halflives_readme.py`: a standard-library Python script that checks the half-lives and the file format;
- `check_E_halflives_readme_6.1.3d.out`: its output on RadioactiveDecay6.1.3d;
- `check_E_halflives_readme.out`: the same on 6.1.2.

## Possible fixes

- Update README_RDM to describe the actual format: three-column summary lines; the five or six columns of the detail lines, including the Q-value in keV; the actual convention of the detail percentages (or renormalise them so that the README is right); the line-length rules and the 120-character limit; the half-life of the `P` lines, ignored for all modes; the list of decay-mode keywords.
- In `LoadDecayTable`, read the lines with `std::getline` into a `std::string`, or at least test `fail()` and issue a G4Exception that names the file, instead of spinning to the loop guard.

## Related reports

Reports on RadioactiveDecay6.1.2 and 6.1.3, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (bug 2780, fixed in 6.1.3)
- B. Ground states with a spurious IT branch, or matched to a floating level (bug 2794, fixed in 6.1.3d)
- D. Decay-mode fractions that ignore the ENSDF branchings
- E. README_RDM does not describe the file format (this report)
- F. Long-lived levels never decayed correctly: no decay data, or E = 0 levels that cannot be created
