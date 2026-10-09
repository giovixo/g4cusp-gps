# RadioactiveDecay6.1.2: inconsistent data make G4RadioactiveDecay emit no beta particles, or no radiation at all, for some nuclides

**Component:** processes/hadronic/models/radioactive_decay, data library RadioactiveDecay6.1.2 (with G4ENSDFSTATE3.0)
**Geant4 version:** 11.4.1 (datasets RadioactiveDecay6.1.2, G4ENSDFSTATE3.0, PhotonEvaporation6.1.2)
**Platform:** macOS 15.8, Apple silicon (arm64)
**Severity:** physics results are wrong, and no warning is printed

## Summary

Some entries of RadioactiveDecay6.1.2 are inconsistent. `G4VRadioactiveDecay::LoadDecayTable` uses these entries as they are. As a result, for some nuclides Geant4 decays the nucleus but emits no beta particle, or emits no radiation at all. No G4Exception or warning is issued in any of these cases. Four classes of problems have a visible effect in Geant4:

1. **Q ≤ 0 in the detail lines (75 levels):** no electron (β⁻) or positron (β⁺) is emitted.
   - Examples: **Mo99**, Ru103, I132, Y92, Y88, Cd105, Cd107, Sn113.
   - In the β⁻ cases, the antineutrino carries the whole decay energy.
2. **Ground states with an IT branch onto themselves (10 levels):** that fraction of decays produces no particle at all.
   - Examples: Y92 50%, Br92 50%, Eu140 50%, Au190 50%, Ir186 20%, Tb154 18%.
3. **Ground states matched to a floating level (6 levels):** the ground state gets the decay modes of a different, floating level, including its IT branch.
   - Example: the Tm164 ground state emits nothing in 80% of its decays.
4. **Levels with no decay data (56 levels that Geant4 can create, T½ ≥ 1 s):** the nucleus is killed when it decays, without emitting anything.
   - Examples: Ta178[0.000X] (2.36 h), Cm252 (2.0 d), Bi198 (10.3 min).

A fifth issue is only in the data files and has no effect in Geant4: 80 half-lives on the `P` lines disagree with ENSDFSTATE (see section 5).

Attached:
- `check_rdm_data.py`: a standard-library Python script that finds all affected levels in the data files;
- `check_rdm_data.out`: its full output;
- `rdecay01_cases.mac`: a macro for the unmodified example `extended/radioactivedecay/rdecay01` that reproduces the Geant4 behaviour.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_rdm_data.py               # data checks; reads $GEANT4_DATA_DIR
cd <rdecay01 build dir>
./rdecay01 rdecay01_cases.mac           # Geant4 behaviour
```

The macro decays 20,000 nuclei of each nuclide at rest, with `/rdecay01/fullChain false` (single decay). It reads the table "Nb of generated particles" of each run.

## Results with rdecay01

Counts of particles generated in 20,000 single decays, from the "Nb of generated particles" table:

| Nuclide | Problem | e⁻ | e⁺ | ν̄ₑ | νₑ | Expected |
|---|---|---|---|---|---|---|
| Y90 (control) | none | 20000 | 0 | 20000 | 0 | as observed |
| Tc99 (control) | none | 20000 | 0 | 20000 | 0 | as observed |
| Co57 (control) | none | 0 | 0 | 0 | 20000 | as observed |
| **Mo99** | 1 | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| Ru103 | 1 | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| I132 | 1 | **0** | 0 | 20000 | 0 | 20000 e⁻ |
| Y88 | 1 | 0 | **0** | 0 | **45** | ~42 e⁺ (0.21%), 20000 νₑ |
| Y92 | 1 + 2 | **0** | 0 | **10007** | 0 | 20000 e⁻, 20000 ν̄ₑ |
| Tb154 | 2 | 0 | 472 | 0 | **16347** | 20000 νₑ |
| Tm164 (ground state) | 3 | 0 | 216 | 0 | **3925** | 20000 νₑ |
| Ta178[0.000X] | 4 | 0 | 0 | 0 | **0** | 20000 νₑ (EC/β⁺ to Hf178) |
| Cm252 | 4 | 0 | 0 | 0 | 0 | 20000 β⁻/α decays |

In every case the decay itself takes place with the correct ENSDFSTATE mean life. The missing products are the β particles, or the whole decay. For issues 1 and 2 the excited daughter levels are still produced, so their de-excitation gammas are still emitted later in a full-chain simulation.

## 1. Q ≤ 0 in the detail lines: no beta particle emitted

In 75 levels, the Q column of non-IT detail lines is zero or negative. It appears to be stored as Q(ground state) − E(daughter level) with Q(ground state) = 0. Example, `z42.a99` (Mo99, real Q_β⁻ = 1357.8 keV):

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

**Effect in Geant4:** for these levels no e⁻ or e⁺ is produced, and in β⁻ decays the antineutrino takes the full Q-value. For Mo99, the antineutrino mean energy is 1.07 MeV with no electron, while the Q-value is 1.36 MeV. The β spectrum, the positron annihilation radiation and the energy deposited by β particles are all missing.

**Affected levels:** Ni59, Rb98, Rb98[270], Y88, Y92, Y94, Y96, Y96[1140], Y97, Y97[667.52], Y97[3522.6], Zr85, Nb84, Nb87, Nb87[3.9], Nb100, Nb100[314], Nb101, Mo87, **Mo99**, Tc86, Tc101, Ru103, Ru105, Ru107, Ru108, Rh98, Rh98[56.3], Rh105, Rh106, Rh106[137], Rh108, Rh111, Rh116, Rh116[150], Pd126, Ag117, Ag117[28.6], Ag118, Ag118[127.63], Ag120, Ag120[203], Cd105, Cd107, Sn106, Sn113, Sn133, Te133, I132, I132[120], Eu130, Tb150, Tb162, Dy163, Ho141[66], Ho154, Ho171, Er165, Yb164, Lu150, Lu156[0X], Hf161, W161, Ir194[190X], Au180, Au181, Bi195, Bi196[271], Po212[2930], Fr214[121], Ra214[1865.2], U228, Pu233, Am234 and Cm248. Most of them lie in the region Z = 37–53. The fraction of decays affected per level is listed in section B of `check_rdm_data.out`.

## 2. Ground states with an IT branch onto themselves

Ten ground states, without a floating flag, have an `IT` summary line. Example, `z39.a92` (Y92, β⁻ 100% in ENSDF):

```
P            0  -        12744
                              BetaMinus            0             0.5
                                     IT            0             0.5
                              BetaMinus            0  -       85.741            0    uniqueFirstForbidden
```

The β⁻ detail lines alone sum to 100%.

**Effect in Geant4:** the IT fraction of decays produces no particle at all. Y92 emits a ν̄ₑ in 10007 of 20000 decays; Tb154 emits a νₑ in 16347 of 20000 decays (expected: every decay).

**Affected levels:**

| Nuclide | IT fraction |
|---|---|
| Br92 | 0.5 |
| Y92 | 0.5 |
| Eu140 | 0.5 |
| Au190 | 0.5 |
| No253 | 0.444 |
| Ir186 | 0.2 |
| Tb154 | 0.179 |
| Sb128 | 0.035 |
| Tc102 | 0.020 |
| Np240 | 0.0012 |

## 3. Ground states matched to a floating level of the library

ENSDFSTATE lists a fixed ground state together with a floating level at the same energy, but the decay library only contains the floating one. `LoadDecayTable` takes "the first level which matches excitation energy regardless of floating level", so the ground state gets the decay table of the other level. Example, Tm164:

```
ENSDFSTATE:   69  164   0  -   1.687953e+11 ns   (T1/2 = 1.95 min, ground state, EC/β+ 100%)
              69  164   0 +X   4.414647e+11 ns   (T1/2 = 5.1 min, isomer)
z69.a164:     P   0 +X   306     ... IT 0.8, EC/β+ 0.2
```

**Effect in Geant4:** the Tm164 ground state emits nothing in 80% of its decays (3925 νₑ in 20000 decays). The decay modes are also those of the wrong level.

**Affected levels:** Rh92, La128, Ho156[52.37], Tm164, Ta160 and Md258. Section D of `check_rdm_data.out` lists the decay modes each one is given.

## 4. Levels without decay data

ENSDFSTATE gives 147 levels with 1 s ≤ T½ < 10¹⁵ s that have no matching level in the decay library. For a floating E = 0 level, or a ground state with no file, `LoadDecayTable` returns an empty decay table. `DecayIt` then kills the nucleus at its decay time without secondaries (the "has no decays" branch). There is no warning, even though the nucleus has a finite ENSDFSTATE lifetime.

I created each of the 147 levels with `/gun/ion Z A 0 E flb` in rdecay01:
- 56 are created as requested and **emit nothing in any decay**;
- 90 cannot be created as requested (Geant4 returns the ground state instead), so they are not affected;
- 1 (Pm152[150X]) decays normally.

The 56 silent levels:

| Level | T½ |
|---|---|
| Bk248[0.000Z] | 9.0 y |
| Cm252 | 2.0 d |
| Lr264 | 4.9 h |
| **Ta178[0.000X]** | 2.36 h |
| Rf267 | 1.3 h |
| Db266 | 24 min |
| Bi198 | 10.3 min |
| Md254 | 10 min |
| Es246[0.000X] | 7.5 min |
| Cm236 | 6.8 min |
| Tl191[0.000X] | 5.2 min |
| Po222 | 2 min |
| Hg187 | 1.9 min |
| 43 more | 1–84 s |

The full list is in the reproducer output. Ta178[0.000X] is produced, for example, in proton-induced spallation of tungsten. In our activation study, 8.5×10⁷ protons (5–700 MeV) on a satellite mass model produced it 224 times, and all of its decay radiation is lost.

**Expected:** decay data for these levels, or at least a warning when an unstable nucleus with a finite ENSDFSTATE lifetime has no decay channel.

## 5. Half-lives on the `P` lines (no effect in Geant4)

The README of RadioactiveDecay says that the half-life on the `P` lines is ignored and that lifetimes are taken from ENSDFSTATE. I confirmed that Geant4 uses the ENSDFSTATE values. However, 80 `P`-line half-lives differ from ENSDFSTATE by more than 1%, and several are exchanged between a ground state and an isomer:

| Nuclide | `P` line | ENSDFSTATE |
|---|---|---|
| Sc44 | 58.6 h | 4.04 h |
| Sc44[271.241] | 4.04 h | 58.6 h |
| Eu152 | 9.31 h | 13.5 y |
| Eu152[45.5998] | 13.5 y | 9.31 h |
| Cu68 / Cu68[721.26] | 3.75 min / 30.9 s | 30.9 s / 3.75 min |
| Se79 | 235 s | 3.27×10⁵ y |
| Tb158[110.3] | 0.4 ms | 10.7 s |

The comment headers have the same problem; for example, `z39.a90` says `# 90Y ( 3.19 H )`. This doesn't affect Geant4, but it misleads anyone who reads the files directly. We used them in an external Bateman-equation code and got wrong activities. It would help to correct these values, or to remove them if they are not meant to be used.

## Possible fixes

- Regenerate the affected entries (issues 1–3) from ENSDF.
- In `LoadDecayTable`:
  - issue a warning, or a JustWarning G4Exception, when the file has no level matching an unstable nucleus;
  - issue one when a level with `noFloat` is matched to a floating level;
  - issue one when a detail line has Q ≤ 0 for a mode other than IT.
