#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"


// Constructor
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    particleGun = new G4GeneralParticleSource();
}


// Destructor
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete particleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    particleGun -> GeneratePrimaryVertex(anEvent);
}


