#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmStandardPhysics.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"

// Constructor
PhysicsList::PhysicsList()
:   fEmPhysicsList(0),
fHadronPhys(),
fDecPhysicsList(0)
{
    SetVerboseLevel(0);
    // Default cut value  (1.0 mm is the Geant4 default)
    defaultCutValue = 100*um;
    
    // Decay physics and all particles
    fDecPhysicsList = new G4DecayPhysics(0);
    fRadioDecPhysicsList = new G4RadioactiveDecayPhysics(0);

    // EM physics (polarization ON)
    // fEmPhysicsList = new G4EmLivermorePhysics(0);
    fEmPhysicsList = new G4EmLivermorePolarizedPhysics(0);

    // Hadronic physics (vector of processes)
    fHadronPhys.push_back( new G4HadronElasticPhysics(0));
    fHadronPhys.push_back( new G4HadronPhysicsFTFP_BERT(0));
//    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BERT(0));
    fHadronPhys.push_back( new G4StoppingPhysics(0));
    fHadronPhys.push_back( new G4IonPhysics(0));
    fHadronPhys.push_back( new G4NeutronTrackingCut(0));
}


// Destructor
PhysicsList::~PhysicsList()
{
    delete fEmPhysicsList;
    for(size_t i=0; i<fHadronPhys.size(); i++) {delete fHadronPhys[i];}
    delete fDecPhysicsList;
}


void PhysicsList::ConstructParticle()
{
    fDecPhysicsList->ConstructParticle();
    fRadioDecPhysicsList->ConstructParticle();
}


void PhysicsList::ConstructProcess()
{
    // In this method, static member functions should be called
    // for all physics processes which you want to use.
    AddTransportation();
    // Electromagnetic physics list
    fEmPhysicsList->ConstructProcess();
    // Hadronic physics lists (activate each element of the vector)
    for(size_t i=0; i<fHadronPhys.size(); i++)
    {
        fHadronPhys[i]->ConstructProcess();
    }
    // Decay physics list
    fDecPhysicsList->ConstructProcess();
    fRadioDecPhysicsList->ConstructProcess();
}



void PhysicsList::SetCuts()
{
    // G4VUserPhysicsListModular::SetCutsWithDefault method sets
    // the default cut value for all particle types
    // (see the class constructor to modify it)
    SetCutsWithDefault();
    // If the verbose level is >1, prints the table of the cut values
    if (verboseLevel>1) DumpCutValuesTable();
}

