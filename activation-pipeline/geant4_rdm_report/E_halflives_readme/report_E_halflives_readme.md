# RadioactiveDecay6.1.3: a 138-character line in `z86.a197` makes `LoadDecayTable` fail (HAD_RDM_100), and README_RDM does not describe the file format

- **Component:** data library RadioactiveDecay6.1.3 (`z86.a197`, README_RDM, comment headers) and processes/hadronic/models/radioactive_decay (`G4VRadioactiveDecay::LoadDecayTable`)
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.3 from bug 2780, final files of 29 October 2026; G4ENSDFSTATE3.0)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** minor. A spurious warning, and one block of data that cannot be read; misleading documentation

This is one of the related reports on the RadioactiveDecay6.1.2 and 6.1.3 data (see "Related reports" at the end).

In 6.1.2, 80 `P`-line half-lives disagreed with ENSDFSTATE. They were corrected in 6.1.3 (bug 2780): the check script finds none. This report covers what is left.

## Summary

1. Line 7 of `z86.a197` (Rn197), a comment added in 6.1.3, is 138 characters long. `LoadDecayTable` reads each line into a 120-character buffer, so the read fails, the loop spins until its 100000-iteration guard and Geant4 issues warning HAD_RDM_100 "While loop count exceeded". The `0 +X` block after that line can never be read.
2. README_RDM describes the summary lines with the wrong number of columns, does not document the Q column of the detail lines (used by Geant4), states a convention for the detail percentages that about half of the levels don't follow, and does not mention the line-length rules, including the 120-character limit of item 1.
3. The half-lives in the comment headers still differ from ENSDFSTATE in 164 files.

Attached:
- `check_E_halflives_readme.py`: a standard-library Python script that checks the half-lives in the files and the file format, and lists the lines of 120 characters or more;
- `check_E_halflives_readme_6.1.3c.out`: its output on the final 6.1.3 files;
- `check_E_halflives_readme.out`: the same on 6.1.2;
- `rdecay01_E_line_length.mac`: a macro for the unmodified example `extended/radioactivedecay/rdecay01` that decays Rn197.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_E_halflives_readme.py <RadioactiveDecay6.1.3 dir>   # default: 6.1.2 from $GEANT4_DATA_DIR
export G4RADIOACTIVEDATA=<RadioactiveDecay6.1.3 dir>
cd <rdecay01 build dir>
./rdecay01 rdecay01_E_line_length.mac
```

## 1. Line of 138 characters in `z86.a197`

```
3:P            0  -        0.065
4:                                  Alpha            0               1
5:                                  Alpha            0  -          100         7411
6:#Modified by hand
7:#Level energy of next parent is 199 in ENSDF decay data  but corresponds to 0+X in ENSDF Levels and gamma  Data, so we take the 0+X energy
8:# to be consistent with PhotoEvaporation and G4ENSDFSTATE data
9:P            0 +X        0.019
10:                                  Alpha            0               1
11:                                  Alpha          100  -          100         7510
```

`LoadDecayTable` reads the file with

```
char inputChars[120]={' '};
...
while (!complete && !DecaySchemeFile.getline(inputChars, 120).eof()) {
```

On a line of 120 characters or more, `getline` stores the first 119 characters and sets `failbit`, not `eofbit`. The loop goes on, but every later `getline` fails at once and returns an empty buffer, so the rest of the file is never read and the loop only ends at its guard of 100000 iterations, with the warning:

```
*** G4Exception : HAD_RDM_100
      issued by : G4VRadioactiveDecay::LoadDecayTable()
While loop count exceeded
```

With rdecay01 (2000 decays of Rn197), the warning is issued when the decay table of Rn197 is loaded. The decay itself is right (α 100%, Q = 7.411 MeV), because the ground-state block comes before line 7. With line 7 shortened to fewer than 120 characters the warning disappears and the "Nb of generated particles" table is identical.

The `0 +X` block (lines 9–11) cannot be read. In Geant4 11.4.1 `/gun/ion 86 197 0 0 X` returns the ground state, because G4IonTable ignores the floating-level base at E = 0 (report F). So this has no visible effect today; if Rn197[0X] could be created, it would get an empty decay table and be killed without emitting anything (report F).

No other line in 6.1.3 has 120 characters or more (the next longest has 115), and none in 6.1.2.

## 2. README_RDM does not match the files

README_RDM is unchanged in 6.1.3.

- It describes the summary lines as four columns (mode, 0, floating flag, fraction). The files have three (mode, 0, fraction): the second column is always 0 and there is no floating flag, neither as a separate column nor attached to another one. The script counts 9845 three-column lines in 6.1.3 (9832 in 6.1.2) and no four-column ones, and `LoadDecayTable` reads exactly three fields (`theDecayMode >> dummy >> decayModeTotal`).
- It describes the detail lines as five columns but names only four fields (mode, daughter level, floating flag, branching ratio). The fifth, the Q-value in keV, is not documented, although `LoadDecayTable` passes it to the decay channel as the Q-value (see bug 2780).
- It says the detail percentages are relative to the mode total. In 1488 of the 3234 levels with detail lines (1486 of 3228 in 6.1.2), the percentages of at least one mode do not add up to 100 (±1); they are relative to all decays. This is harmless while the summary fractions are right, because `LoadDecayTable` rescales the detail lines of each mode to its summary fraction (but see report D).
- `LoadDecayTable` tells summary and detail lines apart by line length (after removing trailing blanks: summary < 72 characters, detail without β type < 84, with β type ≥ 84), not by column count, and it cannot read lines of 120 characters or more (item 1). The README mentions neither rule, so a hand-edited file can easily break them, as `z86.a197` shows.

## 3. Half-lives in the comment headers

The first line of each file gives a half-life, e.g. `# 90Y (  64.05 H  )`. In 6.1.3 those of Y90, Sc44 and Eu152 were corrected, but 164 of the 2195 headers that can be compared with the ENSDFSTATE ground state still differ by more than 1% (194 in 6.1.2). Examples:

| File | Header | ENSDFSTATE (ground state) |
|---|---|---|
| `z101.a258` (Md258) | 57 min | 51.5 d |
| `z98.a237` (Cf237) | 0.8 s | 2.1 s |
| `z98.a240` (Cf240) | 40.3 s | 57.6 s |
| `z100.a245` (Fm245) | 4.2 s | 5.6 s |

Geant4 ignores them, but they mislead anyone who reads the files directly (we used the 6.1.2 values in an external Bateman-equation code and got wrong activities). The full list is in `check_E_halflives_readme_6.1.3c.out`.

## Possible fixes

- Shorten or split line 7 of `z86.a197`, and check in the data production that no line has 120 characters or more.
- In `LoadDecayTable`, read the lines with `std::getline` into a `std::string`, or at least test `fail()` and issue a G4Exception that names the file, instead of spinning to the loop guard.
- Update README_RDM to describe the actual format: three-column summary lines, the Q column, the convention of the detail percentages, and the line-length rules.
- Correct the comment-header half-lives from ENSDFSTATE, or remove them.

## Related reports

Reports on RadioactiveDecay6.1.2 and 6.1.3, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (bug 2780, fixed in 6.1.3)
- B. Ground states with a spurious IT branch, or matched to a floating level: wrong decay modes, some decays emit nothing
- C. Replaced by report F
- D. Wrong decay-mode totals in some levels
- E. A 138-character line in `z86.a197` (6.1.3) breaks `LoadDecayTable`; README_RDM does not match the file format (this report)
- F. Long-lived levels never decayed correctly: no decay data, or E = 0 levels that cannot be created
