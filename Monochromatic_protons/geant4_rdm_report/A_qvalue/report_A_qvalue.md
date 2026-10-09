# RadioactiveDecay6.1.2: Q ≤ 0 in the detail lines of 75 levels, so no β particle is emitted (Mo99, Ru103, I132, Y88, …; regression for Ir194[190X])

- **Component:** processes/hadronic/models/radioactive_decay, data library RadioactiveDecay6.1.2
- **Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2, G4ENSDFSTATE3.0, PhotonEvaporation6.1.2)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** physics results are wrong, and no warning is printed

This is one of five related reports on the RadioactiveDecay6.1.2 data (see "Related reports" at the end).

## Summary

In 75 levels of RadioactiveDecay6.1.2, the Q column of β⁻, β⁺ or EC detail lines is zero or negative. `G4VRadioactiveDecay::LoadDecayTable` passes this column to the decay channel as the Q-value. For these branches Geant4 emits no electron (β⁻) or positron (β⁺), and in β⁻ decay the antineutrino carries all the decay energy. No warning is issued.

The affected levels include **Mo99**, Ru103, I132, Y88, Y92, Cd105, Cd107 and Sn113. For Ir194[190X] this is a **regression** of a fix made in RadioactiveDecay 5.1.1.

Attached:
- `check_A_qvalue.py`: a standard-library Python script that lists all affected levels in the data files;
- `check_A_qvalue.out`: its output;
- `rdecay01_A_qvalue.mac`: a macro for the unmodified example `extended/radioactivedecay/rdecay01` that shows the Geant4 behaviour.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_A_qvalue.py               # reads $GEANT4_DATA_DIR
cd <rdecay01 build dir>
./rdecay01 rdecay01_A_qvalue.mac
```

The macro decays 20,000 nuclei of each nuclide at rest, with `/rdecay01/fullChain false` (single decay). Read the table "Nb of generated particles" of each run.

## Results with rdecay01

Counts of particles generated in 20,000 single decays:

| Nuclide | e⁻ | e⁺ | ν̄ₑ | νₑ | Expected |
|---|---|---|---|---|---|
| Y90 (control) | 20000 | 0 | 20000 | 0 | as observed |
| Tc99 (control) | 20000 | 0 | 20000 | 0 | as observed |
| Co57 (control) | 0 | 0 | 0 | 20000 | as observed |
| **Mo99** | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| Ru103 | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| I132 | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| Y88 | 0 | **0** | 0 | **45** | ~42 e⁺ (0.21%), 20000 νₑ |
| Ir194[190X] | **598** | 0 | 20000 | 0 | 20000 e⁻ |

The decay itself takes place with the correct ENSDFSTATE mean life, and the excited daughter levels are produced, so their de-excitation gammas are still emitted in a full-chain simulation. What is missing is the β particle. For Y88 the EC branches with Q ≤ 0 emit no neutrino either.

For Mo99 the antineutrino mean energy is 1.07 MeV, with no electron, while the decay Q-value is 1.36 MeV. The β spectrum, the positron annihilation radiation and the energy deposited by β particles are all missing.

## The data

The Q column of the detail lines appears to be stored as Q(ground state) − E(daughter level), with Q(ground state) = 0. Example, `z42.a99` (Mo99, real Q_β⁻ = 1357.8 keV):

```
P            0  -     237326.4
                              BetaMinus            0               1
                              BetaMinus     142.6836  -       82.126    -142.6836
                              BetaMinus      509.096  -        1.159     -509.096
                              BetaMinus       920.58  -       16.385      -920.58
```

Correct example, `z39.a90` (Y90):

```
                              BetaMinus            0  -       99.988       2275.6    uniqueFirstForbidden
```

**Affected levels (75):** Ni59, Rb98, Rb98[270], Y88, Y92, Y94, Y96, Y96[1140], Y97, Y97[667.52], Y97[3522.6], Zr85, Nb84, Nb87, Nb87[3.9], Nb100, Nb100[314], Nb101, Mo87, **Mo99**, Tc86, Tc101, Ru103, Ru105, Ru107, Ru108, Rh98, Rh98[56.3], Rh105, Rh106, Rh106[137], Rh108, Rh111, Rh116, Rh116[150], Pd126, Ag117, Ag117[28.6], Ag118, Ag118[127.63], Ag120, Ag120[203], Cd105, Cd107, Sn106, Sn113, Sn133, Te133, I132, I132[120], Eu130, Tb150, Tb162, Dy163, Ho141[66], Ho154, Ho171, Er165, Yb164, Lu150, Lu156[0X], Hf161, W161, Ir194[190X], Au180, Au181, Bi195, Bi196[271], Po212[2930], Fr214[121], Ra214[1865.2], U228, Pu233, Am234 and Cm248.

Most of them lie in the region Z = 37–53. `check_A_qvalue.out` gives, for each level, the number of affected branches, the fraction of decays affected (computed as `LoadDecayTable` does, by rescaling the detail lines of each mode to the mode's summary fraction) and the first Q values.

## Regression for Ir194[190X]

The History file says, for version 5.1.1: "Correction of negative Q value for metastable 190+X level in z77.a194. Bug mentioned by D. Wright." In 6.1.2 the 97% branch of that level again has Q = −20.14 keV:

```
P          190 +X        69048
                              BetaMinus            0               1
                              BetaMinus      2099.55  -            3       318.75
                              BetaMinus      2438.44  -           97       -20.14
```

In rdecay01 only the 3% branch emits an electron (598 of 20,000 decays).

The database was regenerated from ENSDF in version 6.0, so other corrections made in 5.x may have been lost too; version 6.1 already had to restore several files corrected in 5.x.

## Possible fixes

- Regenerate the Q column of the affected levels from ENSDF (Q of the parent level minus the daughter level energy), and check that the corrections made in 5.x are kept.
- In `LoadDecayTable`, issue a warning (or a JustWarning G4Exception) when a detail line has Q ≤ 0 for a mode other than IT.

## Related reports

Five reports on RadioactiveDecay6.1.2, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (this report)
- B. Ground states with a spurious IT branch, or matched to a floating level: wrong decay modes, some decays emit nothing
- C. Unstable levels without decay data are killed without emitting anything, and without a warning
- D. Wrong decay-mode totals in some levels
- E. `P`-line half-lives disagree with ENSDFSTATE; README_RDM does not match the file format

Y92 is affected by both A and B.
