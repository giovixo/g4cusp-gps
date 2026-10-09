# g4cusp-gps

A GEANT4 simulation for CUSP. This branch (`activation`) computes the radioactive nuclides produced in
the CUSP mass model by a primary particle flux.

Source model: General Purpose Source (GPS)

## Requirements

- GEANT4 11.x (tested with 11.4.1), built with GDML support (Xerces-C) and, for the GUI, Qt/OpenGL
- CMake 3.16 or newer
- Python 3 (only for the tools in `tools/`)

## Build

```sh
source <geant4-install>/bin/geant4.sh
cmake -S . -B build
cmake --build build
```

Use `-DWITH_GEANT4_UIVIS=OFF` to build a batch-only executable.
At every build the macros (`macros/`) and the GDML mass model (`gdml-mass-model/`) are copied to the build
directory. Run the program from there: it reads them from the working directory.

## Run

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
without `-n` measures the contribution of neutron capture to the activation.
The two runs write files with the same names, so run them in separate directories.

The primary flux is set with the `/gps/` commands in the macro (see `macros/macrotest.mac`).
The kernel is initialised before the macro runs, so `/run/initialize` in the macro is optional.

User commands:

| Command                        | Default | Description                                              |
|--------------------------------|---------|----------------------------------------------------------|
| `/cusp/stepping/verbose <0/1>` | 1       | print one `*** RADIOISOTOPE` line per recorded nuclide   |

## Output

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

## Documentation

See `doc/` (MkDocs): `cd doc && mkdocs serve`.

## License

MIT, see `LICENSE`.
