# RadioactiveDecay6.1.2: 80 `P`-line half-lives disagree with ENSDFSTATE, and README_RDM does not describe the file format

- **Component:** data library RadioactiveDecay6.1.2 (files `z*.a*` and README_RDM)
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2, G4ENSDFSTATE3.0)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** minor. No effect on Geant4 results; misleading for anyone who reads the files directly

This is one of five related reports on the RadioactiveDecay6.1.2 data (see "Related reports" at the end).

## Summary

1. In 80 levels the half-life on the `P` line differs from ENSDFSTATE by more than 1%. Several are exchanged between a ground state and an isomer (Sc44, Eu152, Cu68). Geant4 ignores these values, but the files present them as data.
2. README_RDM describes the summary lines with the wrong number of columns, does not document the Q column of the detail lines (used by Geant4), and states a convention for the detail percentages that about half of the levels don't follow.

Attached:
- `check_E_halflives_readme.py`: a standard-library Python script that lists the half-life differences and checks the file format;
- `check_E_halflives_readme.out`: its output;
- `rdecay01_E_halflives.mac`: a macro for the unmodified example `extended/radioactivedecay/rdecay01` showing that Geant4 uses the ENSDFSTATE half-lives.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_E_halflives_readme.py     # reads $GEANT4_DATA_DIR
cd <rdecay01 build dir>
./rdecay01 rdecay01_E_halflives.mac
```

## 1. Half-lives on the `P` lines

README_RDM says that the half-life on the `P` lines is ignored for alpha decay, beta decay and IT, and that lifetimes are taken from ENSDFSTATE. In fact `LoadDecayTable` reads it into a dummy variable for every mode, and Geant4 uses the ENSDFSTATE values. However, 80 `P`-line half-lives differ from ENSDFSTATE by more than 1%. Examples:

| Nuclide | `P` line | ENSDFSTATE | Geant4 (rdecay01) |
|---|---|---|---|
| Sc44 | 58.6 h | 4.04 h | 4.04 h (τ = 5.831 h) |
| Sc44[271.241] | 4.04 h | 58.6 h | 58.6 h (τ = 3.523 d) |
| Eu152 | 9.31 h | 13.5 y | 13.5 y (τ = 19.51 y) |
| Eu152[45.5998] | 13.5 y | 9.31 h | 9.31 h (τ = 13.43 h) |
| Cu68 / Cu68[721.26] | 3.75 min / 30.9 s | 30.9 s / 3.75 min | – |
| Se79 | 235 s | 3.27×10⁵ y | 3.27×10⁵ y (τ = 4.721×10⁵ y) |
| Tb158[110.3] | 0.4 ms | 10.7 s | – |

The Geant4 column is the mean life τ printed by rdecay01 for the primary, converted to T½. The full list is in `check_E_halflives_readme.out`.

The comment headers have the same problem; for example, `z39.a90` says `# 90Y ( 3.19 H )` (T½ = 64 h). This doesn't affect Geant4, but it misleads anyone who reads the files directly. We used them in an external Bateman-equation code and got wrong activities. It would help to correct these values, or to remove them if they are not meant to be used.

## 2. README_RDM does not match the files

- It describes the summary lines as four columns (mode, 0, floating flag, fraction). The files have three: there is no floating-flag column. The script counts 9832 three-column lines and no four-column ones.
- It describes the detail lines as five columns but names only four fields (mode, daughter level, floating flag, branching ratio). The fifth, the Q-value in keV, is not documented, although `LoadDecayTable` passes it to the decay channel as the Q-value (see bug 2780).
- It says the detail percentages are relative to the mode total. In 1486 of the 3228 levels with detail lines, the percentages of at least one mode do not add up to 100 (±1); they are relative to all decays. This is harmless while the summary fractions are right, because `LoadDecayTable` rescales the detail lines of each mode to its summary fraction (but see report D).
- `LoadDecayTable` tells summary and detail lines apart by line length (after removing trailing blanks: summary < 72 characters, detail without β type < 84, with β type ≥ 84), not by column count, and reads each line into a 120-character buffer. Every line in 6.1.2 satisfies these rules (the longest is 113 characters), but the README does not mention them, so a hand-edited file can easily break them.

## Possible fixes

- Correct the `P`-line half-lives and comment headers from ENSDFSTATE, or remove them.
- Update README_RDM to describe the actual format: three-column summary lines, the Q column, the convention of the detail percentages, and the line-length rules.

## Related reports

Five reports on RadioactiveDecay6.1.2, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (bug 2780)
- B. Ground states with a spurious IT branch, or matched to a floating level: wrong decay modes, some decays emit nothing
- C. Unstable levels without decay data are killed without emitting anything, and without a warning
- D. Wrong decay-mode totals in some levels
- E. `P`-line half-lives disagree with ENSDFSTATE; README_RDM does not match the file format (this report)
