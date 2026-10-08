# Code

## cmake file

The cmake file `CMakeLists.txt` builds the `cusp-activation` executable.
At every build, the macros (`macros/*.mac`) and the GDML mass model (`gdml-mass-model/*`) are copied
to the build directory, because the program reads them from the working directory.

## c++

Main program: `cusp-activation.cc`

The code works in single or in multi-thread mode. The number of worker threads is set on the command line
(default 1; `0` means all available cores):

```sh
./cusp-activation -t 4 macros/batch.mac
```

If set, the `G4FORCENUMBEROFTHREADS` environment variable overrides `-t`.

The `-n` (`--neutron-activation`) option removes `G4NeutronTrackingCut` from the physics list, so slow
neutrons are tracked until capture. Without it, neutrons are killed 10 µs after they are created.

Header files derived from GEANT4 classes:

* `DetectorConstruction.hh`: reads the GDML mass model and defines the custom materials (e.g. GAGG, FR4)

* `PhysicsList.hh`: Livermore EM, QBBC hadronic physics, decay and radioactive decay; `G4NeutronTrackingCut` unless `-n` is given

* `PrimaryGeneratorAction.hh`: General Particle Source (configured with `/gps/` commands)

* `UserRunAction.hh`, `UserRun.hh`, `UserEventAction.hh`

* `SteppingAction.hh`: records the radioactive nuclides produced in the geometry

For example `DetectorConstruction` inherits from `G4VUserDetectorConstruction`

## Output

Each nuclide with a lifetime between 0.1 s and 1e18 s is recorded when it is produced, and then killed
(it is not decayed). The ntuple `Events` is written as CSV, one file per run and per worker thread:
`scorefile_run<N>_nt_Events_t<thread>.csv`. A thread that records no nuclide writes no file.

| Column     | Type   | Description                                  |
|------------|--------|----------------------------------------------|
| `EventID`  | int    | Event number                                 |
| `Isotope`  | string | Nuclide name (e.g. `Al26`)                   |
| `Lifetime` | double | Mean lifetime in seconds                     |
| `Volume`   | string | Physical volume where the nuclide was produced |

To merge the thread files of each run into `scorefile_run<N>.csv` (sorted by `EventID`):

```sh
python3 tools/merge_scorefiles.py <output-dir> [--run N] [-o OUTDIR] [--delete]
```

The command `/cusp/stepping/verbose 0` switches off the `*** RADIOISOTOPE` line printed for each
recorded nuclide (default 1).
