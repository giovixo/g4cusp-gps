# g4cusp-gps

A GEANT4 simulation for CUSP. This branch (`activation`) computes the activation background of the CUSP
CubeSat polarimeter: the radioactive nuclides produced in the CUSP mass model by trapped protons, their
activities along the orbit, and the rates and spectra their decays give in the detector. It follows the
method of Campana et al. 2026 (Exp. Astron. 61:21).

The repository has two GEANT4 programs and a Python pipeline:

| Part | What it does |
|---|---|
| `cusp-activation` | protons on the mass model → production yields of the nuclides, per volume |
| `cusp-postactivation` | decays one nuclide in one volume → energy deposits in the 96 scintillators |
| `Monochromatic_protons/` | Python pipeline (steps 0–8): runs both programs, solves the decay chains, folds with the orbit flux |

Source model: General Purpose Source (GPS)

## Requirements

- GEANT4 11.x (tested with 11.4.1), built with GDML support (Xerces-C) and, for the GUI, Qt/OpenGL
- CMake 3.16 or newer
- Python ≥ 3.11 with numpy, scipy, pandas, matplotlib, mpmath (for the pipeline and the tools in `tools/`)

## Build

```sh
source <geant4-install>/bin/geant4.sh
cmake -S . -B build
cmake --build build
```

This builds `cusp-activation` and `cusp-postactivation`, which share the geometry and the physics list
(static library `cusp-common`). Use `-DWITH_GEANT4_UIVIS=OFF` to build batch-only executables.
At every build the macros (`macros/`) and the GDML mass model (`gdml-mass-model/`) are copied to the build
directory. Run the programs from there: they read them from the working directory.

Sources: `src/` and `include/`, each split into `common/` (geometry, physics list), `activation/` and
`postactivation/` (classes prefixed `PostAct`).

Physics: Livermore EM, `G4HadronElasticPhysicsXS`, QBBC inelastic, stopping, `G4IonPhysicsXS`, decay and
radioactive decay.

## `cusp-activation`

```sh
cd build
./cusp-activation                        # interactive GUI
./cusp-activation macrotest.mac          # batch mode
./cusp-activation -t 4 macrotest.mac     # batch mode, 4 worker threads (0 = all cores)
./cusp-activation -n macrotest.mac       # batch mode, neutron activation on
```

By default `G4NeutronTrackingCut` kills neutrons 10 µs after they are created, before most of them
thermalise. With `-n` (`--neutron-activation`) the cut is removed, so slow neutrons are tracked until
they are captured or escape. Everything else in the physics list stays the same, so running with and
without `-n` measures the contribution of neutron capture to the activation (less than 3e-5 of the yields
for CUSP). The two runs write files with the same names, so run them in separate directories.

The primary flux is set with the `/gps/` commands in the macro (see `macros/macrotest.mac`).
The source must lie inside the world volume (a 50 cm box), or the primaries cross it without interacting.
The standard source is a sphere of radius 13 cm centred at (0, −1.5, 4.8) cm, cosine law, θmax = 90°.
The kernel is initialised before the macro runs, so `/run/initialize` in the macro is optional.

User commands:

| Command                        | Default | Description                                              |
|--------------------------------|---------|----------------------------------------------------------|
| `/cusp/stepping/verbose <0/1>` | 1       | print one `*** RADIOISOTOPE` line per recorded nuclide   |

### Output

Each nuclide with a lifetime between 0.1 s and 1e18 s is recorded where it is produced and then killed
(it is not decayed). The ntuple `Events` is written as CSV, one file per run and per worker thread:

```
scorefile_run<N>_nt_Events_t<thread>.csv
```

Columns: `EventID` (int), `Isotope` (e.g. `Al26`), `Lifetime` (mean lifetime, s), `Volume` (physical volume).

Merge the thread files into one file per run (`scorefile_run<N>.csv`, sorted by `EventID`):

```sh
python3 ../tools/merge_scorefiles.py .            # all runs in the current directory
python3 ../tools/merge_scorefiles.py . --run 0    # only run 0
```

The pipeline (`0_run.py`) does this merge itself.

## `cusp-postactivation`

```sh
./cusp-postactivation [-t nthreads] [macro]
```

One run decays one nuclide, at rest and uniformly distributed in one volume, and records the energy
deposited in the scintillators. It is normally driven by `Monochromatic_protons/run_postactivation.py`
(step 4), which writes the macros.

| Command                        | Default | Description                                                    |
|--------------------------------|---------|----------------------------------------------------------------|
| `/postact/isotope <name>`      |         | nuclide, with the pipeline names (e.g. `Na22`, `Co60[58.590]`, `Ta178[0.000X]`) |
| `/postact/volume <name>`       |         | physical volume where the nuclide decays                       |
| `/postact/output <prefix>`     | postact | prefix of the output files                                     |
| `/postact/minHalfLife <t> <unit>` | 1 us | daughters with T½ ≥ t are killed before they decay             |
| `/postact/timeWindow <t> <unit>`  | 10 us | deposits later than t after the decay are dropped            |

Time 0 is the decay. Output:

- `<prefix>_t<N>.csv`: `RunID,EventID,ScintID,Edep_keV,t_ns`. ScintID 0–63 = `PV-Scatterer_001`–`064`
  (plastic), 64–95 = `PV-Absorber_001`–`032` (GAGG).
- `<prefix>_runs.csv`: `RunID,Volume,Isotope,NDecays,G4Ion`.

The first run of a process with a given prefix deletes the old files of that prefix; later runs append.

Note: GEANT4 creates only one ground-state level per nuclide, so a few isomers named with an energy of 0
(e.g. `Ta178`, `Pm136`) decay as their `[0.000X]` level. The `G4Ion` column records the ion that was used.

## Activation pipeline

`Monochromatic_protons/` holds the pipeline. `Monochromatic_protons/PIPELINE_GUIDE.md` documents every
step and has a complete example run.

| Step | Script | Output |
|---|---|---|
| 0 | `0_run.py` + `config_cusp.toml` | GEANT4 activation runs, one per proton energy |
| 1 | `activation_parser.py` | `results.pkl` (yields per primary) |
| 2 | `build_decay_chains.py` | `DecayChains/` |
| 3 | `compute_activities.py` | `activities.pkl`, `active_isotopes.pkl` |
| 4 | `run_postactivation.py` | `result_postact/` (hits per decay) |
| 4b | `build_spectra.py` | `result_postact/spectra_<mode>.npz` (event selection) |
| 5 | `spenvis_parser.py` | orbit flux plots (optional) |
| 6 | `accumulate_spectra.py` | `spectra.pkl`, `count_rate.dat` |
| 7 | `activation_history.py` | count-rate history along the orbit |
| 8 | `average_spectrum.py` | mean out-of-belt spectrum with line identification |
| report | `build_report.py` | HTML and PDF report (one or more orbits) |

Steps 0–4 do not depend on the orbit: a new orbit or flux model (SPENVIS file) needs only steps 6–8, or
`build_report.py`.

The decay data are read from the GEANT4 data libraries: source `geant4.sh` before steps 2–3 (the tools use
`$G4RADIOACTIVEDATA` or `$GEANT4_DATA_DIR`).

`Monochromatic_protons/geant4_rdm_report/` holds reports on errors found in the GEANT4 radioactive decay
data (report A is GEANT4 bug 2780).

## Documentation

- `Monochromatic_protons/PIPELINE_GUIDE.md`: the activation pipeline.
- `doc/` (MkDocs): `cd doc && mkdocs serve`.

## License

MIT, see `LICENSE`.
