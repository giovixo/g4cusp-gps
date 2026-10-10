#ifndef PostActDetectorConstruction_h
#define PostActDetectorConstruction_h 1

#include "DetectorConstruction.hh"

// Mass model of the common code plus the sensitive detector on the 96 scintillators.
// ScintID: Scatterer_NNN -> NNN-1 (0-63), Absorber_NNN -> 64+NNN-1 (64-95).

class PostActDetectorConstruction : public DetectorConstruction
{
public:
    PostActDetectorConstruction();
    virtual ~PostActDetectorConstruction();

    // Called once per worker thread (and in the sequential mode)
    virtual void ConstructSDandField();

    // Work-around for a thread-safety bug of Geant4 11.4.1: G4Voxelizer::GetCandidates reads its
    // mutable std::map with operator[], which inserts the missing entries, so two threads calling
    // Inside() on the same G4TessellatedSolid can corrupt it (crash "pointer being freed was not
    // allocated" in G4TessellatedSolid::InsideVoxels). Inserting every entry once, on the master
    // after the geometry is built and before the workers start, makes the later calls read-only.
    // Returns the number of tessellated solids processed.
    static G4int PrefillVoxelCandidates();
};

#endif
