#ifndef PostActSteppingAction_h
#define PostActSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class PostActEventAction;

// 1. Time origin: when the primary decays, the global time of its products is shifted
//    by -t0 (t0 = time of the decay), so that time 0 is the decay. The primary can decay
//    after 1e18 s, and a double could not resolve nanoseconds on such a time.
// 2. Daughters with half-life >= PostActConfig::MinHalfLife() are killed when they stop
//    (their recoil energy deposit is kept): their activity is computed separately.

class PostActSteppingAction : public G4UserSteppingAction
{
public:
    explicit PostActSteppingAction(PostActEventAction* eventAction);
    virtual ~PostActSteppingAction();

    virtual void UserSteppingAction(const G4Step* step);

private:
    PostActEventAction* fEventAction;
};

#endif
