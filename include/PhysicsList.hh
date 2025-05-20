#ifndef PhysicsList_H
#define PhysicsList_H 1

#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

class PhysicsList: public G4VModularPhysicsList
{
public:
    PhysicsList();        // Constructor
    ~PhysicsList();        // Destructor
    
protected:
    // Construct particle and physics processes
    void ConstructParticle();
    void ConstructProcess();
    
    // Define cuts
    void SetCuts();
    
private:
    G4VPhysicsConstructor*               fEmPhysicsList;
    std::vector<G4VPhysicsConstructor*>  fHadronPhys;
    G4VPhysicsConstructor*               fDecPhysicsList;
    G4VPhysicsConstructor*               fRadioDecPhysicsList;

};

#endif








