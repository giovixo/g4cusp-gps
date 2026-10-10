# CLAUDE.md

Geant4 simulation of the CUSP CubeSat polarimeter, branch `activation`: it computes the radioactive nuclides
produced in the CUSP mass model by trapped protons, decays them in the detector and turns them into an
activation background (rates and spectra). It follows the method of Campana et al. 2026 (Exp. Astron. 61:21).
`Campana_2026.pdf` in the root is git-ignored on purpose (publisher copyright): never commit it.

## Environment and build

- Geant4 11.4.1 at `~/geant4-v11.4.1`. Run `source ~/geant4-v11.4.1/bin/geant4.sh` before building, running
  Geant4 or running the Python steps that need decay data (they find it through `$GEANT4_DATA_DIR`).
- The repo root holds the user's **in-source Xcode build** (`CMakeCache.txt`, `cusp-activation.xcodeproj`,
  `Debug/`, root copies of `*.mac`/`*.gdml`/`*.xml`). These are not stale: leave them alone. One Xcode project,
  two schemes; `Debug/cusp-activation` and `Debug/cusp-postactivation` are the defaults of `0_run.py` and
  `run_postactivation.py` (`--executable` overrides; the step-4 production used a Release scratch build).
- For a scratch build use `cmake -S . -B <dir>`. A bare `cmake <srcdir>` reconfigures the in-source cache.
  Keep `WITH_GEANT4_UIVIS=ON`.

## Branches

- **Never touch `main`**: no commits, merges, rebases, pushes or checkouts, and no worktrees or branches based on
  it. The Agent tool's `isolation: "worktree"` bases worktrees on `main`: do not use it here; for parallel agents,
  create worktrees yourself from `activation`, or give agents disjoint files in the main tree.
- **Work on `activation` only**: it holds the whole pipeline, steps 0–8 (the post-activation work was merged into
  it by fast-forward on 2026-10-10). `post-activation` is retired: the local branch is deleted; the remote one
  (created by the repository owner) is kept with one extra commit, `BRANCH_RETIRED.md`, a reminder that it can be
  deleted (`git push origin --delete post-activation`) once the owner agrees.

## Geant4 applications

- Sources: `src/` + `include/` split into `common/` (`DetectorConstruction`, `PhysicsList`, built as the static
  library `cusp-common`), `activation/` and `postactivation/` (classes prefixed `PostAct`).
- Physics: Livermore EM, `G4HadronElasticPhysicsXS`, QBBC inelastic, stopping, `G4IonPhysicsXS`, decay +
  radioactive decay (no `_HP` neutron models).
- Geometry: GDML in `gdml-mass-model/`, 173 tessellated volumes, all direct daughters of the world, no rotations.
  World is a 50 cm box; the payload bounding box is x −49…47, y −70…40, z −46…142 mm. Detector: 64 plastic
  scatterers `PV-Scatterer_NNN` and 32 GAGG absorbers `PV-Absorber_NNN`.

### `cusp-activation` (steps 0–1)

- `cusp-activation [-t nthreads] [-n] [macro]`: `-t 0` = all cores; `-n` removes `G4NeutronTrackingCut`
  (neutron activation). No macro opens the GUI.
- `SteppingAction` records each nuclide with mean life 0.1 s–1e18 s where it is produced and kills it: the output
  is production yields only (`EventID, Isotope, Lifetime, Volume`, one CSV per thread).
- **The source must lie inside the world volume**, or primaries cross it without interacting. Standard source:
  sphere R = 13 cm centred at (0, −1.5, 4.8) cm, cosine law, θmax = 90° (beam area πR²sin²θmax = 530.9 cm²).
  See `macros/macrotest.mac`.

### `cusp-postactivation` (step 4)

- `cusp-postactivation [-t nthreads] [macro]`. Commands: `/postact/isotope` (step-3 names, e.g.
  `Ta178[0.000X]`), `/postact/volume`, `/postact/output`, `/postact/minHalfLife` (1 us, = step 2
  `min_halflife`), `/postact/timeWindow` (10 us). One run = one (volume, isotope) pair.
- Ion at rest, uniform in the volume (global transform, rejection sampling, no fallback). The ion is created on
  the master when `/postact/isotope` is given; level energies are matched within the 1 eV rounding of the
  pipeline names. Time 0 = the decay; daughters with T½ ≥ minHalfLife are killed when they stop.
- The RDM very-long-decay threshold is set to 1e60 y in `main` (the Geant4 default, 1 y, would stop Na22 decays).
- Output: `<prefix>_t<N>.csv` (`RunID,EventID,ScintID,Edep_keV,t_ns`; ScintID 0–63 = PV-Scatterer_001–064,
  64–95 = PV-Absorber_001–032, matched by exact name and material) and `<prefix>_runs.csv`
  (`RunID,Volume,Isotope,NDecays,G4Ion`). The first run of a process with a prefix deletes that prefix's old
  files; runs within a process append.

### Geant4 behaviour to know

- Geant4 creates only **one E = 0 level per nuclide** (`G4IonTable::FindIon` falls back to any E = 0 ion):
  `Ta178`, `Ta178[0.000Y]`, `Pr134`, `Pm135`, `Pm136`, `Re172` decay as their `[0.000X]` level (warning
  PostAct022; the `G4Ion` column records it). Step 1 only ever produced the `[0.000X]` levels; the other names
  come from the chain builder, which therefore does not fully mimic Geant4 for these daughters.
- The chain builder deliberately reproduces Geant4's decay-data behaviour, data errors included (see
  `activation-pipeline/geant4_rdm_report/`; report A is Geant4 bug 2780).
- One step-4 production batch (1 in 67) aborted once with a heap error reported in
  `G4TessellatedSolid::InsideVoxels`; it did not reproduce (rerun, AddressSanitizer, ThreadSanitizer, 1.7e8
  concurrent `Inside()` calls). It is **not** a demonstrated G4Voxelizer bug. The driver reruns incomplete
  batches, so outputs stay complete. A work-around based on a wrong diagnosis (`3d7c1c0`) was reverted (`0c532b8`).

## Activation pipeline (`activation-pipeline/`)

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
| 8 | `average_spectrum.py` | out-of-belt mean spectrum |
| report | `build_report.py` + `report_template.html` | `report_run/cusp_activation_report.{html,pdf}` |

- Step 0 production (17 energies × 1e8 protons, all cores) takes about 4 h; run it in the background under
  `caffeinate -i`. `0_run.py` skips energies already done with the same nprim unless `--overwrite` is given.
- Neutron activation (`-n`) changes the yields by less than 3e-5 at all energies: `output/` without `-n` is the
  reference set.
- Photonuclear activation (`G4EmExtraPhysics`: gamma-, electro-, muon-nuclear), tested 2026-10-10 with 1e7 protons
  at 20, 60, 150, 300, 700 MeV and the creator process of every nuclide: 11 photonNuclear nuclides out of 2.5e5
  (1.2e-4 of the yield at 700 MeV, from pi0 photons; 0 at 150 and 300 MeV). Weighted with the orbit spectra:
  ~5e-5 of the production (AP8MIN 4e-5, AP9 6e-5; < 3e-4 at 95% CL). Not included in the physics list.
- Step 4 production (2026-10-10): 3e8 decays, N ∝ out-of-belt activity (pilot: N ∝ A·√p within 10%), 22,787
  pairs, about 2 h on 11 cores; the 2472 pairs below the step-3 threshold (~0.9% of the activity) were
  deliberately not simulated. `run_postactivation.py` resumes: complete batches are skipped, incomplete ones rerun.
- Step 4b: the event selection is **provisional** (thresholds 5 keV plastic / 20 keV GAGG, 500 ns window, no
  resolution or quenching; modes `scat_single`, `abs_single`, `compton` = exactly 1 + 1, `any`). Rebuilding the
  spectra takes ~20 s. Steps 6 and 8 read them with `--spectra-file result_postact/spectra_<mode>.npz` and
  `--emin/--emax`; the old `result_spectra/*.dat` input still works.
- The user is mostly interested in **out-of-belt** results: averages in steps 7–8 and in the report use
  out-of-belt steps by default (`--belt-threshold`, `--all-orbit`).
- SPENVIS files (`AP[89]*.txt`) are git-ignored: `AP8MIN.AP8.output_mean_flux_550km_SSO.txt` (stops at 400 MeV)
  and `AP9MEAN.AP9.output_mean_flux_550km_SSO.txt` (to 2000 MeV, uses the 700 MeV runs).
- **A new orbit or flux model needs only steps 6–8** (or `build_report.py` with several SPENVIS files and
  `--labels`, one report with ratios): steps 0–4 do not depend on the orbit, so never rerun them for that.
- Steps 4 (weights), 7 and 8 use F(>E_min), E_min = lowest simulation energy (recorded by step 6 in
  `count_rate.dat`/`spectra.pkl`), as the flux time profile; the belt passages use the total flux.
- `myUtilities` (plot style) is on the user's own PYTHONPATH, not in the repo; its import is optional.

## Conventions

- Commit or push only when asked. Remote `https://github.com/giovixo/g4cusp-gps.git` belongs to another user;
  Riccardo pushes as a collaborator over HTTPS (keychain token). Commits are signed with an SSH key from
  1Password: if `git commit` fails with a 1Password error, ask the user to unlock it; never disable signing.
- Keep, unless asked otherwise: `SteppingAction` of `cusp-activation` (no stacking action), the empty `UserRun`
  and `UserEventAction`, the unused materials in `DetectorConstruction`, the fixed random seed.
- Simulation outputs, pickles, logs, plots, reports and `.dat` files stay out of git (see the `.gitignore` files).
- `HANDOFF.md` (git-ignored, local only) holds the current session state and next steps; it is loaded below
  when present.

@HANDOFF.md
