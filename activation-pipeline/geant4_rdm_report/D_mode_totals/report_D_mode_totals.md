# RadioactiveDecay6.1.2 and 6.1.3: decay-mode fractions that ignore the ENSDF branchings (Pu233 decays by α only, Np231 α 95%, Bi212[239] has no α branch, Tb157 feeds the wrong level)

- **Component:** data library RadioactiveDecay6.1.2 and 6.1.3 (read by `G4VRadioactiveDecay::LoadDecayTable`)
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2, G4ENSDFSTATE3.0, PhotonEvaporation6.1.2)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** physics results are wrong, and no warning is printed

This is one of the related reports on the RadioactiveDecay6.1.2 and 6.1.3 data (see "Related reports" at the end). The problem is the same in RadioactiveDecay6.1.3d (provided after bug 2794): the mode fractions of the levels below are unchanged.

## Summary

Each level in a RadioactiveDecay file has one summary line per decay mode (the fraction of decays in that mode) and detail lines (the daughter levels of each mode, with intensities in %). `LoadDecayTable` rescales the detail lines of each mode so that they add up to the summary fraction. **Only the summary lines therefore decide how often each mode occurs.**

In ENSDF, the decay-mode branchings of a level (for example %ε and %α) and the feeding of individual daughter levels are separate data. The branchings are totals, often measured directly. The feeding of a daughter level is known only where the decays to that level have been observed, through their γ rays or α lines. For many nuclides, especially short-lived ones, only part of a mode's decays, or none of them, can be assigned to definite daughter levels in this way. The detail lines of the RadioactiveDecay files contain only these assigned feedings, so for such a mode they add up to less than its branching.

In some levels the summary fractions are not the ENSDF branchings. They are the sum of the detail intensities of each mode, divided by the sum over all modes. So a mode whose decays are only partly assigned to daughter levels gets too few decays, a mode with no assigned decays disappears, and the other modes are inflated. For example:

- **Pu233** (ENSDF: ε 99.88%, α 0.12%). None of the ε decays is assigned to a daughter level, so the file has no ε line at all, and the α branch, with 0.12% in its detail line, becomes 100% of the decays.
- **Np231** (ENSDF: ε ≈ 98%, α ≈ 2%). Detail lines: α 2%, β⁺ 0.109%. Summary: α 2 / 2.109 = 0.948, β⁺ 0.052.

A script finds 21 candidate levels with T½ ≥ 1 s. I compared each of them with the ENSDF adopted levels in NuDat 3 (see the table): 17 have mode fractions that disagree with ENSDF (two of them, Ho152 and No255, by 2–3σ), one feeds the wrong daughter level (Tb157), and three are consistent with ENSDF or have no ENSDF branching to compare with. Five of the candidates have Z > 100; by default Geant4 does not decay them (`G4NucleusLimits`: Z ≤ 100, A ≤ 260), so they matter only if the limits are raised.

## Examples from the data (6.1.3d)

```
z94.a233 (Pu233; ENSDF: ε 99.88%, α 0.12%)
P            0  -         1254
                                  Alpha            0               1
                                  Alpha           15  -         0.12         6405

z83.a212 (Bi212[239]; ENSDF: α 67%, β⁻ 33%)
P          239  -         1500
                              BetaMinus            0               1
                              BetaMinus            0  -           33       2490.5
P        250.7  -            0
                                  Alpha            0            0.67
                              BetaMinus            0            0.33
                                  ...
```

The modes of Bi212[239] are correct in a second record, at 250.7 keV with T½ = 0. `LoadDecayTable` matches the ENSDFSTATE level at 239 keV to the first record, so this level always decays by β⁻.

**Tb157** is a different case: the mode is right (EC 100%), but the ground-state feeding has intensity 1e-25, so after rescaling every decay goes to Gd157[54.536]. In ENSDF this level gets only 0.11% of the decays (Iγ(54.5 keV) = 0.0084%):

```
z65.a157 (Tb157)
P            0  -  2.24059e+09
                               MshellEC            0               1
                               MshellEC            0  -        1e-25           60
                               MshellEC       54.536  -         0.11        5.464
```

Only M-shell capture is listed; capture to the ground state (Q = 60 keV) is also possible from the K and L shells.

## The 21 candidates compared with ENSDF

Geant4 = mode fractions in RadioactiveDecay6.1.2 and 6.1.3d (identical). ENSDF = adopted levels in NuDat 3 (accessed 4 November 2026). ε stands for EC + β⁺.

| Level | T½ | Geant4 | ENSDF (NuDat 3) | Result |
|---|---|---|---|---|
| Rh98[56.3] | 3.6 min | ε 100% | IT 89%, ε 11% | wrong: no IT |
| Tb157 | 71 y | EC 100%, all to Gd157[54.5] | ε 100%, 99.9% to the ground state | wrong daughter level |
| Ho152 | 2.70 min | α 18.6%, ε 81.4% | α 12 (3)%, ε 88 (3)% | α 2σ too high |
| Tm155[41] | 45 s | α 100% | ε > 98%, α < 2% | wrong |
| Os181[49.2] | 2.7 min | β⁺ 75%, IT 25% | ε ≈ 100%, IT ≤ 3% | wrong: IT too high |
| Tl181 | 3.2 s | α 100% | α ≤ 10%, rest ε | wrong |
| Tl183 | 6.9 s | α 100% | ε > 0, α not known | wrong: no ε |
| Pb185 | 6.3 s | α 100% | α 34 (25)%, ε not known | wrong: no ε |
| Bi212[239] | 25 min | β⁻ 100% | α 67%, β⁻ 33% | wrong: no α |
| Bi212[1478] | 7 min | β⁻ 100% | IT > 75%, β⁻ < 25% | wrong: no IT |
| Np231 | 48.8 min | α 94.8%, β⁺ 5.2% | ε ≈ 98%, α ≈ 2% | wrong |
| Pu233 | 20.9 min | α 100% | ε 99.88%, α 0.12% | wrong |
| Am232 | 1.32 min | α 100% | ε ≈ 97%, α ≈ 3% | wrong |
| Am233 | 3.2 min | α 100% | ε 95.5%, α 4.5% | wrong |
| Cm234 | 51 s | α 100% | α ≈ 74%, ε ≈ 22%, SF ≈ 4% | wrong: no ε, no SF |
| Bk234 | 2.33 min | β⁺ 100% | no branching given (identified by its α decay) | cannot be checked |
| Md248 * | 7 s | α 32%, β⁺ 68% | α 20 (10)%, ε 80 (10)% | consistent within 1.2σ |
| No255 * | 3.52 min | α 44%, β⁺ 56% | α 30 (5)%, ε 70 (5)% | α 2.7σ too high |
| Lr262 * | 4 h | SF 100% | SF < 10%, ε and α not known | wrong |
| Db262 * | 35 s | α 100% | α ≈ 67%, SF + ε ≈ 33% | wrong |
| Sg265[152X] * | 8.5 s | α 100% | α ≥ 49%, SF ≤ 51% | consistent |

\* Z > 100: not decayed with the default nucleus limits.

When the ε branch is kept, it is sometimes all β⁺, with no EC lines (Os181[49.2], Np231, Md248, No255).

## Results with rdecay01

Counts of particles generated in 20,000 single decays (`rdecay01_D_mode_totals.mac`):

| Nuclide | Observed, 6.1.2 | Observed, 6.1.3d | Expected (ENSDF) |
|---|---|---|---|
| Co57 (control) | 20000 νₑ | 20000 νₑ | as observed |
| Pu233 | 20000 U229[15], no α¹ | 20000 α (U229[15]) | about 20000 νₑ (EC), 24 α |
| Np231 | 18922 α, 1078 e⁺ | 18922 α, 1078 e⁺ | about 19600 νₑ, 400 α |
| Am233 | 20000 α | 20000 α | about 19100 νₑ, 900 α |
| Tl183 | 20000 α | 20000 α | mostly νₑ (EC/β⁺) |
| Bi212[239] | 20000 e⁻ (Po212) | 20000 e⁻ (Po212) | about 13400 α, 6600 e⁻ |
| Rh98[56.3] | 15331 νₑ, no e⁺¹ | 15331 e⁺, 20000 νₑ | about 17800 IT, 2200 νₑ |
| Tb157 | 20000 Gd157[54.5] | 20000 Gd157[54.5] | about 20000 Gd157 ground state, 22 Gd157[54.5] |

¹ The difference between 6.1.2 and 6.1.3d comes from the Q values fixed for bug 2780, not from the mode fractions.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_D_mode_totals.py [RadioactiveDecay dir]   # default: RadioactiveDecay6.1.2 in $GEANT4_DATA_DIR
cd <rdecay01 build dir>
./rdecay01 rdecay01_D_mode_totals.mac                   # set G4RADIOACTIVEDATA to choose the library
```

The script lists the levels (1 s ≤ T½ < 10¹⁵ s) where a mode gets more than 50% of the decays although its detail lines add up to less than 50%, or where a mode with 1% or more in its detail lines has no summary line. It is a heuristic: it does not read ENSDF branchings, so levels with similar errors can be missed.

Attached:
- `check_D_mode_totals.py`: standard-library Python script that lists the candidates;
- `check_D_mode_totals.out`, `check_D_mode_totals_6.1.3d.out`: its output for 6.1.2 and 6.1.3d (same 21 levels);
- `rdecay01_D_mode_totals.mac`: macro for the unmodified example `extended/radioactivedecay/rdecay01`.

## Possible fixes

- Take the summary fractions from the ENSDF decay-mode branchings (the `%` values of the adopted level), not from the sums of the daughter-level feedings.
- When none of a mode's decays is assigned to a daughter level (for example the ε branch of Pu233), keep its summary line, with the decay going to the daughter ground state.
- In the data generator, check that each summary fraction agrees with the ENSDF branching within its uncertainty, and that the ground-state feeding is not replaced by a placeholder intensity (Tb157).

## Related reports

Reports on RadioactiveDecay6.1.2 and 6.1.3, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (bug 2780, fixed in 6.1.3)
- B. Ground states with a spurious IT branch, or matched to a floating level (bug 2794, fixed in 6.1.3d)
- D. Decay-mode fractions that ignore the ENSDF branchings (this report)
- E. A 138-character line in `z86.a197` (6.1.3) breaks `LoadDecayTable`; README_RDM does not match the file format (the line and the header half-lives are fixed in 6.1.3d; README_RDM is unchanged)
- F. Long-lived levels never decayed correctly: no decay data, or E = 0 levels that cannot be created
