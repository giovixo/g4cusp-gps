# Code

## cmake file 
The cmake file `CMakeLists.txt` has the following option to activate the debug messages:
```cmake
#----------------------------------------------------------------------------
# Enable the DEBUG directive for this target
#
target_compile_definitions(cusp_gdml PRIVATE DEBUG)
```
Comment to to switch off the debug.

The debug option can be implemented as in the following code example:
```cpp
    #ifdef DEBUG
    std::ostringstream strValue;
    strValue << "Event ID: " << event.event_id << ", "
        << "Time: " << event.time << ", "
        << "Angle: " << event.angle << ", "
        << "Pol_x: " << event.pol_x << ", "
        << "Pol_y: " << event.pol_y << ", "
        << "Pol_z: " << event.pol_z << ", "
        << "Energy: " << event.energy;
    testOutput.print(strValue.str());
    #endif
```

## c++

Main program: `cusp_gdml.cc`

The code works in single or in multi-thread mode. The number of threads can be defined programmatically.

```cpp
// Construct the run manager
#ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(6);
#else
    G4RunManager * runManager = new G4RunManager;
#endif
```

Header files derived from GEANT4 classes:

* `DetectorConstruction.hh`

* `DetectorHit.hh`

* `PhysicsList.hh`

* `PrimaryGeneratorAction.hh`

* `SensitiveDetector.hh`

* ... 

For example `DetectorConstruction`  inherits from `G4VUserDetectorConstruction`

Other header files:

* `WriteToFile.hh`

* `RsmSource.hh`

## How to test the user c++ code outsite the geant4 application

Edit `CMakeLists.txt` in `cpp-lab/<_dir_name>` if necessary

Compile:

```c
cd cpp-lab/<_dir_name>/build
cmake ..
make
```

Run:

`./<program_name>`