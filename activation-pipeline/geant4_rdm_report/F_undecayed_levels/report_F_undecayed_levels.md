# Radioactive decay: 154 long-lived levels are never decayed correctly, and no warning is given. 55 have no decay data and emit nothing; 99 can never be created, because G4IonTable ignores the floating-level base at E = 0 (W178 → Ta178 emits nothing after the EC)

- **Components:** processes/hadronic/models/radioactive_decay (`G4VRadioactiveDecay::LoadDecayTable`, `DecayIt`, data RadioactiveDecay6.1.2 and 6.1.3) and particles (`G4IonTable::GetIon`, `G4IonTable::FindIon`, `G4NuclideTable::GetIsotope`)
- **Geant4 version:** 11.4.1 (datasets G4ENSDFSTATE3.0, RadioactiveDecay6.1.2 and 6.1.3, PhotonEvaporation6.1.2)
- **Platform:** macOS 15.8, Apple silicon (arm64)
- **Severity:** decay radiation is lost, or comes from another level with another lifetime; no warning

This is one of the related reports on the RadioactiveDecay6.1.2 and 6.1.3 data (see "Related reports" at the end). It replaces report C ("Unstable levels without decay data are killed without emitting anything"), announced in the list of report B: part of what C described as harmless has the same cause as a defect in G4IonTable, so the two are reported together.

The numbers are for the final RadioactiveDecay6.1.3 files (29 October 2026, bug 2780); those for 6.1.2 are in `check_F_undecayed_levels.out`. The Geant4 behaviour is the same with both.

## Summary

For Geant4 to decay a level correctly, two lookups must succeed:

1. **Identity:** `G4IonTable::GetIon(Z, A, E, flb)` must create the requested level.
2. **Data:** `G4VRadioactiveDecay::LoadDecayTable` must find the decay data of that level in RadioactiveDecay. Energy and floating-level base must both match.

Each lookup fails, without a warning, for a set of long-lived levels (1 s ≤ T½ < 10¹⁵ s):

- **No decay data (55 levels).** 145 levels have no decay data. 55 of them can be created, and they are killed at their decay time without any secondary: Ta178[0X] (2.36 h), Cm252 (2.0 d), Bi198 (10.3 min) and others.
- **Cannot be created (99 levels).** Some nuclides have more than one level at E = 0: a ground state plus an isomer of unknown excitation energy (`0+X`), or two states of unknown order (`0+X`, `0+Y`). At E = 0, G4IonTable ignores the floating-level base and always returns the first E = 0 level in `ENSDFSTATE.dat`. 172 E = 0 levels in 161 nuclides can never be created; 99 of them have T½ ≥ 1 s. A request for one of them gets the first level instead, with that level's lifetime and decay data. The other 89 levels without data are among these 99, so they are affected through the substitution. Examples: Sb120 (8-, 5.76 d) becomes the 15.9 min ground state; Np240 (7.22 min) becomes the 1.03 h ground state.

The two failures can combine. W178 (T½ 21.6 d) decays 100 % by EC to Ta178[0Y] (1+, 9.31 min), which has decay data but cannot be created. Geant4 creates Ta178[0X] (7-, 2.36 h) instead, which has no decay data. **All Ta178 radiation of the W178 chain is lost**: Hf X-rays, the 93 keV gamma and the β⁺ branch. Bk248 is a second case: the 23.7 h level has data, but becomes the 9 y level, which has none.

Three paths in Geant4 request levels that cannot be created:
- **radioactive decay:** daughters are created with the floating-level base of the decay record (e.g. `G4ECDecay.cc:63`). Eight decays of parents with T½ ≥ 1 s feed such a level;
- **hadronic de-excitation:** `G4PhotonEvaporation` keeps the floating index of the final level (`G4PhotonEvaporation.cc:585-586`), and `G4ExcitationHandler` passes it to `GetIon` (`G4ExcitationHandler.cc:554-561`). Spallation and capture residues that end on a shadowed isomer are therefore created as the first level;
- **user requests,** e.g. `/gun/ion 73 178 0 0 Y`.

Attached:
- `check_F_undecayed_levels.py`: a standard-library Python script. It lists (1) the levels without decay data, (2) the E = 0 levels that cannot be created, and (3) the decay records that feed them, with the consequence of each. It can also write an rdecay01 macro for the levels of list 1 and classify the rdecay01 output;
- `check_F_undecayed_levels_6.1.3c.out`: its output on the final 6.1.3 files, including the classification of the rdecay01 run;
- `check_F_undecayed_levels.out`: the same on 6.1.2;
- `rdecay01_F_undecayed_levels.mac`: a short macro for the unmodified example `extended/radioactivedecay/rdecay01`;
- `rdecay01_F_all_levels.mac`: the macro for the 145 levels of list 1 (6.1.3), written by the script.

## Floating levels in ENSDF

A level energy is measured through the transitions that link it to the ground state. When no transition links a group of levels to the rest of the scheme, ENSDF gives their energies relative to an unknown base: `E+X`, then Y, Z, U, V, W… for further independent groups (ENSDF manual, field E). G4ENSDFSTATE3.0 and RadioactiveDecay keep only the number as the energy and give the letter in a separate column. So `Ta178 0+X` and `Ta178 0+Y` both have E = 0 and are told apart only by the letter. NUBASE2020 estimates the 1+ level of Ta178 at 100 ± 50 keV (from systematics) and flags the order as uncertain.

## How to reproduce

```
source <geant4-install>/bin/geant4.sh
python3 check_F_undecayed_levels.py <RadioactiveDecay6.1.3 dir>    # default: 6.1.2 from $GEANT4_DATA_DIR
export G4RADIOACTIVEDATA=<RadioactiveDecay6.1.3 dir>
cd <rdecay01 build dir>
./rdecay01 rdecay01_F_undecayed_levels.mac                          # the cases below

python3 check_F_undecayed_levels.py <RadioactiveDecay6.1.3 dir> --macro=rdecay01_F_all_levels.mac
./rdecay01 rdecay01_F_all_levels.mac > all_levels.log
python3 check_F_undecayed_levels.py <RadioactiveDecay6.1.3 dir> --log=all_levels.log
```

The short macro decays 2000 nuclei of each nuclide at rest (single decays, then two full chains); the long one decays 200 nuclei of each level. For each run, read the name and mean life of the primary ion and the table "Nb of generated particles".

## Results with rdecay01 (short macro; identical with 6.1.2 and 6.1.3)

| Run | Requested | Geant4 primary | Result | Expected |
|---|---|---|---|---|
| Control | Co57 | Co57 | 2000 νₑ, Fe57 levels | as observed |
| No data | Ta178[0X] (7-, 2.36 h) | `Ta178[0.00000X]`, mean life 3.405 h | **no secondaries** | EC/β⁺ to Hf178 |
| No data | Cm252 (2.0 d) | Cm252 | **no secondaries** | α/β⁻ decays |
| No data | Bi198 (10.3 min) | Bi198 | **no secondaries** | EC/β⁺ to Pb198 |
| Shadowed | `/gun/ion 73 178 0 0 Y` (1+, 9.31 min) | `Ta178[0.00000X]`, mean life 3.405 h | **no secondaries** | EC/β⁺ to Hf178, data in `z73.a178` (`P 0 +Y`) |
| Shadowed | `/gun/ion 51 120 0 0 X` (8-, 5.76 d) | `Sb120`, mean life 22.92 min | decays as the 1+ ground state | the 8- isomer |
| Full chain | W178 | `W178 ---> Ta178[0.00000X]` | **the chain stops**: only the 2000 νₑ of the W178 EC | W178 → Ta178 (1+) → Hf178 |
| Full chain | U240 | `U240 ---> Np240 ---> …` | all 2000 decays through the Np240 ground state (mean life 1.488 h) | 72 % through Np240 (7.22 min) |

## 1. Levels without decay data

The script finds 145 levels with 1 s ≤ T½ < 10¹⁵ s and no matching level in the library (147 in 6.1.2). For a ground state with no file, or a floating E = 0 level that is not in the file, `LoadDecayTable` returns an empty decay table. `DecayIt` then kills the nucleus at its decay time without secondaries (the "has no decays" branch). There is no warning, although the nucleus has a finite ENSDFSTATE lifetime.

Each of the 145 levels was created with `/gun/ion Z A 0 E flb` in rdecay01:
- 55 are created as requested and **emit nothing in any decay** (56 in 6.1.2);
- 89 cannot be created as requested: they are all shadowed E = 0 levels (section 2);
- 1 (Pm152[150X]) decays normally.

The 55 silent levels, by half-life:

| Level | T½ |
|---|---|
| Bk248[0Z] | 9.0 y |
| Cm252 | 2.0 d |
| Lr264 | 4.9 h |
| **Ta178[0X]** | 2.36 h |
| Rf267 | 1.3 h |
| Db266 | 24 min |
| Bi198 | 10.3 min |
| Md254 | 10 min |
| Es246[0X] | 7.5 min |
| Cm236 | 6.8 min |
| Tl191[0X] | 5.2 min |
| Po222 | 2 min |
| 43 more | 1–84 s |

In 6.1.3, Hg187 and Ho154[0X] gained decay data; the others are unchanged. The Ta178[0X] data existed in RadioactiveDecay5.6 (`P 0 +X 8496`, EC to the 1147.4 and 1479.0 keV levels of Hf178) and are missing from 6.1.2 and 6.1.3.

## 2. E = 0 levels that cannot be created

From `ENSDFSTATE.dat`, restricted to the levels that G4NuclideTable loads:
- 161 nuclides have more than one E = 0 level. For 136 of them the first is the fixed ground state, for 25 a floating level (Ta178, Pr134, Re172 and Pm136 have no fixed ground state at all).
- 172 E = 0 levels are shadowed, i.e. never created. 99 of them have T½ ≥ 1 s, and 10 of those 99 have their own decay data in the library, which is never used.

Shadowed levels with T½ ≥ 10 min (all 99 are in the script output):

| Shadowed level | T½ | Own data | Geant4 creates | T½ | Its data |
|---|---|---|---|---|---|
| Sb120[0X] | 5.76 d | no | Sb120 | 15.9 min | yes |
| Bk248[0X] | 23.7 h | **yes** | Bk248[0Z] | 9 y | **no** |
| Tb154[0Y] | 22.7 h | no | Tb154 | 21.5 h | yes |
| Re182[0X] | 14.1 h | no | Re182 | 2.68 d | yes |
| Tb154[0X] | 9.4 h | no | Tb154 | 21.5 h | yes |
| Es256[0X] | 7.6 h | no | Es256 | 25.4 min | yes |
| Es250[0X] | 2.22 h | no | Es250 | 8.6 h | yes |
| Ir186[0X] | 1.9 h | no | Ir186 | 16.6 h | yes |
| Md258[0X] | 57 min | yes | Md258 | 51.5 d | yes |
| Hg191[0X] | 50.8 min | no | Hg191 | 49 min | yes |
| Bi200[0X] | 31 min | no | Bi200 | 36.4 min | yes |
| Md254[0X] | 28 min | no | Md254 | 10 min | no |
| Am246[0X] | 25 min | no | Am246 | 39 min | yes |
| Bi198[0X] | 11.6 min | no | Bi198 | 10.3 min | no |
| Pr134[0Y] | 11 min | no | Pr134[0X] | 17 min | yes |
| **Ta178[0Y]** | 9.31 min | **yes** | Ta178[0X] | 2.36 h | **no** |

Decay records of parents with T½ ≥ 1 s whose daughter is a shadowed level:

| Parent | T½ | Requested daughter | Geant4 daughter | BR | Consequence |
|---|---|---|---|---|---|
| W178 | 21.6 d | Ta178[0Y], 9.31 min | Ta178[0X], 2.36 h | 100 % | no radiation |
| U240 | 14.1 h | Np240[0X], 7.22 min | Np240, 1.03 h | 72 % | other lifetime and decay data |
| Pu246 | 10.8 d | Am246[0X], 25 min | Am246, 39 min | 2 % | other lifetime and decay data |
| At220 | 3.7 min | Bi216[0X], 6.6 min | Bi216, 2.25 min | 8 % | other lifetime and decay data |
| Tm154[70] | 8.1 s | Ho150[0X], 23.5 s | Ho150, 1.2 min | 99.8 % | other lifetime and decay data |
| Zr102 | 2.9 s | Nb102[0X], 1.3 s | Nb102, 4.3 s | 58.5 % | other lifetime and decay data |
| Nd154 | 25.9 s | Pm154[0X], 1.73 min | Pm154, 2.68 min | 5.4 % | other lifetime and decay data |
| Ir173 | 9 s | Re169[0X], 15.1 s | Re169, 8.1 s | 4.1 % | no data for either level |

The hadronic path cannot be counted from the data files. It applies to every shadowed level on which photon evaporation can end, e.g. Sb120 (8-) produced on Sn and Sb, or the Re182 and Ir186 isomers produced in spallation of heavy targets.

## Cause

**No decay data.** For a nucleus with a finite ENSDFSTATE lifetime and no matching `P` record, `LoadDecayTable` returns an empty decay table, and `DecayIt` kills the nucleus without secondaries. Neither issues a warning.

**Cannot be created.** At E = 0 the floating-level base is dropped at three points:
1. `G4NuclideTable::GetIsotope` (`G4NuclideTable.cc:118`) accepts a level if `flb == level->GetFloatLevelBase() || E == 0.0`. At E = 0 it therefore returns the first level of the multimap, i.e. the first in `ENSDFSTATE.dat`.
2. `G4IonTable::GetIon` (`G4IonTable.cc:483-490`): when no ion exists yet, it takes the floating base of that level (`flb1`) and creates the ion with it (`CreateIon(Z, A, E, flb1)`, line 500).
3. `G4IonTable::FindIon` (`G4IonTable.cc:644-648`): "rerun search on ground level without check on floating". At E = 0 it returns the first existing ion with E = 0, whatever its floating base.

All E = 0 levels also get the same PDG encoding. `G4NuclideTable` raises the isomer level only for E > 0 (`if (ionE > 0) ++iLevel`, `G4NuclideTable.cc:211`), and `CreateIon` sets `lvl = 0` when `Eex == 0` (`G4IonTable.cc:277`), so the encoding computed at line 287 is the same for all of them.

`LoadDecayTable`, by contrast, requires the floating base of the `P` record to match that of the ion (`G4VRadioactiveDecay.cc:611-618`). So the data written for a shadowed level are never used, and the level substituted for it may have no data at all.

## Effect in an application

In our activation study of the CUSP CubeSat (5–700 MeV protons on a satellite mass model, Geant4 11.4.1), W178 is produced in the tungsten collimators. Ta178 (1+) is in secular equilibrium with it: 0.345 Bq out of 35.7 Bq of payload activity after one year in orbit. All of its radiation is lost. We reran these decays with a private copy of `z73.a178` in which the `P 0 +Y` block is relabelled `+X` (`/process/had/rdm/setRadioactiveDecayFile`). The missing radiation is 2.6 % of the Compton-event background and 1.6 % of the plastic-scintillator background. Ta178[0X] (0.044 Bq) cannot be corrected this way, because no data exist for it.

## Possible fixes

- In `G4NuclideTable::GetIsotope` and `G4IonTable::FindIon`, ignore the floating base at E = 0 only when the request has none (`flb == no_Float`). A request with a letter should match the level with that letter.
- Give each E = 0 floating level its own isomer level (and PDG encoding), as for E > 0 levels.
- Add decay data for the levels of section 1. For Ta178[0X], the RadioactiveDecay5.6 block could be restored.
- Issue a warning, or a JustWarning G4Exception once per nuclide:
  - in `LoadDecayTable` or `DecayIt`, when an unstable nucleus with a finite ENSDFSTATE lifetime has no decay channel;
  - in `GetIon`, when the level returned has a floating base other than the one requested.

## Related reports

Reports on RadioactiveDecay6.1.2 and 6.1.3, submitted separately:
- A. Q ≤ 0 in the detail lines: no β particle emitted (bug 2780, fixed in 6.1.3)
- B. Ground states with a spurious IT branch, or matched to a floating level: wrong decay modes, some decays emit nothing
- C. Replaced by this report
- D. Wrong decay-mode totals in some levels
- E. A 138-character line in `z86.a197` (6.1.3) breaks `LoadDecayTable`; README_RDM does not match the file format
- F. Long-lived levels never decayed correctly: no decay data, or E = 0 levels that cannot be created (this report)
