#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "WriteToFile.hh"


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

    // Retrieve the primary vertex
    #ifdef DEBUG
    G4int eventID = anEvent->GetEventID();
    G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex();
    // Retrieve the primary particle
    G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
    // Retrieve the energy
    G4double energy = primaryParticle->GetKineticEnergy()/CLHEP::keV;
    // Print the energy
    //testOutput.print("The energy is:");
    testOutput.print(std::to_string(eventID));
    testOutput.print(std::to_string(energy));
    #endif
}


