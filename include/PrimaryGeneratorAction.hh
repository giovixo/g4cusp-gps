
#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

#include "WriteToFile.hh"

// #include "RsmSource.hh"
#include "GlobalRsmSource.hh"

// Mandatory user class that defines the properties of the
// particle gun and run initialization
// Derived from the G4VUserPrimaryGeneratorAction initialisation
// abstract base class.

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();		// Constructor
    virtual ~PrimaryGeneratorAction();					// Destructor
    
public:
    // This method generates the primary particles
    void GeneratePrimaries(G4Event* anEvent);
    
    
private:
    G4ParticleGun* particleGun;
    //RsmSource rsmSource; // to read the source data
    void GenerateVertex(G4double source_h, G4double source_r, G4double source_alpha, G4ThreeVector &vertex_center, G4ThreeVector& vertex); // Generate the source vertex
    void GenerateDirection(G4double detector_h, G4ThreeVector vertex, G4ThreeVector& direction); // Generate the source direction
    void GenerateSquare(G4double number);     // Just a dummy method
    WriteToFile testOutput; // to write the test output
};

#endif


