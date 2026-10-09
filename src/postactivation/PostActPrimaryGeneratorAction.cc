#include "PostActPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"


PostActPrimaryGeneratorAction::PostActPrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);
    fParticleGun->SetParticleDefinition(G4Geantino::Definition());
    fParticleGun->SetParticleEnergy(1*MeV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}


PostActPrimaryGeneratorAction::~PostActPrimaryGeneratorAction()
{
    delete fParticleGun;
}


void PostActPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
