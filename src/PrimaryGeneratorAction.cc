#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"


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


