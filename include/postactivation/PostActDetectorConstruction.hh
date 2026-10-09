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
};

#endif
