#ifndef PostActStackingAction_h
#define PostActStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

class G4ParticleDefinition;

// Kills the long-lived nuclides produced with exactly zero kinetic energy.
// Geant4 runs the at-rest decay of such a track in its first step, before the stepping
// action is called, so PostActSteppingAction could not stop it.

class PostActStackingAction : public G4UserStackingAction
{
public:
    PostActStackingAction();
    virtual ~PostActStackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track);

    // Unstable nucleus with half-life >= PostActConfig::MinHalfLife(): its decay belongs to
    // another run (its activity is computed separately), so it must not decay in the event
    static G4bool IsLongLived(const G4ParticleDefinition* particle);
};

#endif
