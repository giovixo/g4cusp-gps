# CLAUDE.md

Geant4 simulation of the CUSP CubeSat polarimeter, branch `activation`: it computes
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

- **Never touch `main`**: no commits, merges, rebases, pushes or checkouts, and no worktrees or branches
  based on it. The Agent tool's `isolation: "worktree"` bases worktrees on `main`: do not use it here.
- **Work on `activation` only**: it holds the whole pipeline, steps 0–8 (the post-activation work was merged
  into it by fast-forward on 2026-10-10). `post-activation` is retired: the local branch is deleted; the remote
  one (created by the repository owner) is kept with one extra commit, `BRANCH_RETIRED.md`, a reminder that it can
  be deleted (`git push origin --delete post-activation`) once the owner agrees.

## Geant4 applications

- Two executables from one CMake project (one Xcode project, two schemes):
  `cusp-activation` (steps 0–3: nuclide production) and `cusp-postactivation` (step 4: decay of each
  nuclide in each volume, detector response).
- `cusp-postactivation [-t nthreads] [macro]`, commands `/postact/isotope` (step-3 names, e.g. `Ta178[0.000X]`),
  `/postact/volume`, `/postact/output`, `/postact/minHalfLife` (1 us, = step 2 `min_halflife`),
  `/postact/timeWindow` (10 us). Ion at rest, uniform in the volume; time 0 = the decay; daughters with
  T½ ≥ minHalfLife are killed when they stop. Output: `<prefix>_t<N>.csv` (`RunID,EventID,ScintID,Edep_keV,t_ns`,
  ScintID 0–63 = PV-Scatterer_001–064, 64–95 = PV-Absorber_001–032) and `<prefix>_runs.csv`. The first run of a
  process with a prefix deletes that prefix's old files, so a rerun replaces them; runs within a process append.
- Geant4 creates only **one E = 0 level per nuclide** (`G4IonTable::FindIon` falls back to any E = 0 ion):
  `Ta178`, `Ta178[0.000Y]` and `Pr134` decay as `Ta178[0.000X]`/`Pr134[0.000X]` (warning PostAct022; the runs
  file has a `G4Ion` column with the level that decayed). Step 1 only ever produced the `[0.000X]` levels; the
  other names come from the chain builder, which therefore does not fully mimic Geant4 for these daughters.
  The RDM very-long-decay threshold is set to 1e60 y in `main` (default 1 y would stop Na22 decays).
- One production batch (1 in 67) aborted once with a heap error reported in `G4TessellatedSolid::InsideVoxels`;
  it did not reproduce (rerun, AddressSanitizer, ThreadSanitizer, 1.7e8 concurrent `Inside()` calls). It is
  **not** a demonstrated G4Voxelizer bug: `Inside()` does not modify the voxel candidate map. The driver reruns
  incomplete batches, so outputs stay complete. The pre-fill work-around based on that wrong diagnosis
  (`3d7c1c0`) was reverted in `0c532b8`.
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
| 4 | `run_postactivation.py` + `cusp-postactivation` | `result_postact/` (hit lists per decay) |
| 4b | `build_spectra.py` (event selection) + `pair_spectra.py` (format) | `result_postact/spectra_<mode>.npz` |
| 5 | `spenvis_parser.py` | orbit flux plots (optional) |
| 6 | `accumulate_spectra.py` | `spectra.pkl`, `count_rate.dat` |
| 7 | `activation_history.py` | rate history, out-of-belt running average |
| 8 | `average_spectrum.py` | steady-state out-of-belt spectrum |
| report | `build_report.py` + `report_template.html` | `report_run/cusp_activation_report.{html,pdf}` |

- A production set (17 energies × 1e8 protons, all cores) takes about 4 h; run it in the background
  under `caffeinate -i`. `0_run.py` skips energies already done with the same nprim unless `--overwrite` is given.
- Step 4 production (2026-10-10): 3e8 decays, N ∝ out-of-belt activity, 22,787 pairs; the 2472 pairs below the
  step-3 threshold (~0.9% of the activity) were deliberately not simulated. Event selection in `build_spectra.py`
  is provisional (thr 5 keV plastic / 20 keV GAGG, 500 ns window, no resolution or quenching; modes
  `scat_single`, `abs_single`, `compton` = exactly 1+1, `any`). Steps 6 and 8 read it with
  `--spectra-file result_postact/spectra_<mode>.npz` and `--emin/--emax`.
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
