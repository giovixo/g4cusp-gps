#ifndef PostActPrimaryGeneratorAction_h
#define PostActPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

// Primary generator of cusp-postactivation.
// Placeholder: shoots a geantino from the origin. It will generate an ion at rest,
// uniformly distributed in a given physical volume.

class PostActPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PostActPrimaryGeneratorAction();
    virtual ~PostActPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* anEvent);

private:
    G4ParticleGun* fParticleGun;
};

#endif
