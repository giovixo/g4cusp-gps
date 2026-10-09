# CLAUDE.md

Geant4 simulation of the CUSP CubeSat polarimeter, branches `activation` and `post-activation`: it computes
the radioactive nuclides produced in the CUSP mass model by trapped protons and turns them into an
activation background spectrum. It follows the method of Campana et al. 2026 (Exp. Astron. 61:21).
`Campana_2026.pdf` in the root is git-ignored on purpose (publisher copyright): never commit it.

## Environment and build

- Geant4 11.4.1 at `~/geant4-v11.4.1`. Run `source ~/geant4-v11.4.1/bin/geant4.sh` before building,
  running Geant4 or running pipeline steps 2–3 (they find the decay data through `$GEANT4_DATA_DIR`).
- The repo root holds the user's **in-source Xcode build** (`CMakeCache.txt`, `cusp-activation.xcodeproj`,
  `Debug/`, root copies of `*.mac`/`*.gdml`/`*.xml`). These are not stale: leave them alone.
  `Debug/cusp-activation` is the executable used by `Monochromatic_protons/0_run.py`.
- For a scratch build use `cmake -S . -B <dir>`. A bare `cmake <srcdir>` reconfigures the in-source cache.
  Keep `WITH_GEANT4_UIVIS=ON`.

## Branches

- **Never touch `main`** (no commits, merges or pushes).
- `post-activation` = `activation` (merged in `9896885`) + the post-activation program (step 4). Work happens
  here; when step 4 works, `activation` is fast-forwarded to it and the two branches become one.

## Geant4 applications

- Two executables from one CMake project (one Xcode project, two schemes):
  `cusp-activation` (steps 0–3: nuclide production) and `cusp-postactivation` (step 4: decay of each
  nuclide in each volume, detector response; work in progress).
- Sources are split into `src/` + `include/` subdirectories: `common/` (`DetectorConstruction`, `PhysicsList`,
  built as the static library `cusp-common`), `activation/` and `postactivation/` (classes prefixed `PostAct`).
- `cusp-activation [-t nthreads] [-n] [macro]`: `-t 0` = all cores; `-n` removes `G4NeutronTrackingCut`
  (neutron activation). No macro opens the GUI.
- Physics: Livermore EM, `G4HadronElasticPhysicsXS`, QBBC inelastic, stopping, `G4IonPhysicsXS`,
  decay + radioactive decay (no `_HP` neutron models).
- `SteppingAction` records each nuclide with mean life 0.1 s–1e18 s where it is produced and kills it:
  the output is production yields only (`EventID, Isotope, Lifetime, Volume`, one CSV per thread).
- Geometry: GDML in `gdml-mass-model/`. World is a 50 cm box; the payload bounding box is
  x −49…47, y −70…40, z −46…142 mm.
- **The source must lie inside the world volume**, or primaries cross it without interacting.
  Standard source: sphere R = 13 cm centred at (0, −1.5, 4.8) cm, cosine law, θmax = 90°
  (beam area πR²sin²θmax = 530.9 cm²). See `macros/macrotest.mac`.

## Activation pipeline (`Monochromatic_protons/`)

`PIPELINE_GUIDE.md` documents every step and has a complete example run.

| Step | Script | Output |
|---|---|---|
| 0 | `0_run.py` + `config_cusp.toml` | `output/scorefile_<E>_MeV_<nprim>.csv`, `run_info.json` |
| 1 | `activation_parser.py` | `results.pkl` |
| 2 | `build_decay_chains.py` | `DecayChains/` |
| 3 | `compute_activities.py` | `output/activities.pkl`, `active_isotopes.pkl`, … |
| 4 | Geant4 post-activation run (another code, not in this repo) | `result_spectra/` |
| 5 | `spenvis_parser.py` | orbit flux plots (optional) |
| 6 | `accumulate_spectra.py` | `spectra.pkl`, `count_rate.dat` |
| 7 | `activation_history.py` | rate history, out-of-belt running average |
| 8 | `average_spectrum.py` | steady-state out-of-belt spectrum |

- A production set (17 energies × 1e8 protons, all cores) takes about 4 h; run it in the background
  under `caffeinate -i`. `0_run.py` skips energies already done with the same nprim unless `--overwrite` is given.
- The user is mostly interested in **out-of-belt** results: averages in steps 7–8 use out-of-belt steps
  by default (`--belt-threshold`, `--all-orbit`).
- SPENVIS files (`AP[89]*.txt`) are git-ignored. The current one,
  `AP8MIN.AP8.output_mean_flux_550km_SSO.txt`, stops at 400 MeV; the 700 MeV runs are kept on purpose
  for AP9 spectra.
- Neutron activation (`-n`) changes the yields by less than 3e-5 at all energies (production run of
  2026-10-09): `output/` without `-n` is the reference set.
- The chain builder deliberately reproduces Geant4's decay-data behaviour, data errors included
  (see `geant4_rdm_report/bug_report.md`).
- `myUtilities` (plot style) is on the user's own PYTHONPATH, not in the repo; its import is optional.

## Conventions

- Commit or push only when asked. Remote `https://github.com/giovixo/g4cusp-gps.git` belongs to another
  user; Riccardo pushes as a collaborator over HTTPS (keychain token).
- Keep, unless asked otherwise: `SteppingAction` (no stacking action), the empty `UserRun` and
  `UserEventAction`, the unused materials in `DetectorConstruction`, the fixed random seed.
- Simulation outputs, pickles, logs, plots and `.dat` files stay out of git (see the `.gitignore` files).
- `HANDOFF.md` (git-ignored, local only) holds the current session state and next steps; it is loaded
  below when present.

@HANDOFF.md
