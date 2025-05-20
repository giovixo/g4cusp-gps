#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "WriteToFile.hh"

// Mandatory user class that defines the properties of the
// particle gun and run initialization
// Derived from the G4VUserPrimaryGeneratorAction initialisation
// abstract base class.

class G4GeneralParticleSource;
class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();		// Constructor
    virtual ~PrimaryGeneratorAction();					// Destructor
    
public:
    // This method generates the primary particles
    void GeneratePrimaries(G4Event* anEvent);
    
    
private:
    G4GeneralParticleSource*    particleGun;
    WriteToFile testOutput; // to write the test output
};

#endif


