<!-- DRAFT: check the ENSDF branching ratios quoted below against NuDat before submitting. -->

# RadioactiveDecay6.1.2: wrong decay-mode totals in some levels (Pu233 emits no α and never decays by EC, Np231 α 95%, Tl183 α 100%, Tb157 always feeds Gd157[54.5])

- **Component:** data library RadioactiveDecay6.1.2 (read by `G4VRadioactiveDecay::LoadDecayTable`)
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2, G4ENSDFSTATE3.0, PhotonEvaporation6.1.2)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** physics results are wrong, and no warning is printed

This is one of five related reports on the RadioactiveDecay6.1.2 data (see "Related reports" at the end).

## Summary

README_RDM says that the detail percentages are relative to the total of their mode. `LoadDecayTable` follows this: it rescales the detail lines of each mode so that they add up to the mode's summary fraction. Only the summary lines therefore decide how often each mode occurs.

In some levels the summary fractions are wrong: they give a mode far more decays than ENSDF, or leave a mode out. They seem to come from renormalising the detail percentages to 1, so that a mode with no placed intensity disappears and the other modes are inflated. A heuristic finds 21 candidate levels with T½ ≥ 1 s. Four of them were checked with rdecay01: Pu233, Np231, Tl183 and Tb157.

Attached:
- `check_D_mode_totals.py`: a standard-library Python script that lists the candidate levels;
- `check_D_mode_totals.out`: its output;
- `rdecay01_D_mode_totals.mac`: a macro for the unmodified example `extended/radioactivedecay/rdecay01` that shows the Geant4 behaviour.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_D_mode_totals.py          # reads $GEANT4_DATA_DIR
cd <rdecay01 build dir>
./rdecay01 rdecay01_D_mode_totals.mac
```

The macro decays 20,000 nuclei of each nuclide at rest, with `/rdecay01/fullChain false` (single decay). Read the table "Nb of generated particles" of each run.

## Results with rdecay01

20,000 single decays each:

| Nuclide | Observed | Expected (ENSDF) |
|---|---|---|
| Co57 (control) | 20000 νₑ, Fe57 excited levels | as observed |
| Pu233 | U229[15] in every decay, **no α particle**¹ | EC 99.88% to Np233, α 0.12% |
| Np231 | α 18922 (95%), β⁺ 1078 (5%) | EC/β⁺ 98%, α 2% |
| Tl183 | α 20000 (100%) | mostly EC/β⁺, α a few % |
| Tb157 | every decay is M-shell EC to **Gd157[54.536]** | EC 100%, to the ground state |

¹ The α detail line has Q = −15 keV, so no α particle is emitted (see report A).

## The data

Two examples:

```
z94.a233 (Pu233; ENSDF: EC 99.88%, α 0.12%)
P            0  -         1254
                                  Alpha            0               1
                                  Alpha           15  -         0.12          -15

z81.a183 (Tl183)
                                  Alpha            0               1
                               BetaPlus            0               0
   detail lines: Alpha 50 %, BetaPlus 50 %
```

In Pu233 the EC branch is missing and the α branch, with 0.12% of decays in its detail line, becomes 100% of the decays. Its Q column is also −15 keV, so the α particle is not emitted (see report A). In Tl183 the β⁺ branch has detail lines but a summary fraction of 0.

**Candidate levels (21):** Rh98[56.3], Tb157, Ho152, Tm155[41], Os181[49.2], Tl181, Tl183, Pb185, Bi212[239], Bi212[1478], Np231, Pu233, Am232, Am233, Cm234, Bk234, Md248, No255, Lr262, Db262 and Sg265[152X].

They were selected by a heuristic: the summary gives more than 50% of the decays to a mode whose detail lines add up to less than 50%, or about 0 to a mode whose detail lines add up to 1% or more. `check_D_mode_totals.out` lists the summary fractions and the detail sums of each one. Some of them may be correct, so each case should be checked against ENSDF.

## Possible fixes

- Regenerate the summary lines of these levels from the ENSDF decay-mode branchings, not from the detail intensities.
- In the data generator, check that each summary fraction agrees with the corresponding ENSDF branching.

## Related reports

Five reports on RadioactiveDecay6.1.2, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted
- B. Ground states with a spurious IT branch, or matched to a floating level: wrong decay modes, some decays emit nothing
- C. Unstable levels without decay data are killed without emitting anything, and without a warning
- D. Wrong decay-mode totals in some levels (this report)
- E. `P`-line half-lives disagree with ENSDFSTATE; README_RDM does not match the file format
