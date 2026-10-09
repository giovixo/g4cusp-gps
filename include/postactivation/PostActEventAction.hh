#ifndef PostActEventAction_h
#define PostActEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class PostActRunAction;
class PostActScintillatorSD;

// At the end of each event writes one row per scintillator with a non-zero deposit.
// The totals are read from the (thread-local) sensitive detector.

class PostActEventAction : public G4UserEventAction
{
public:
    explicit PostActEventAction(PostActRunAction* runAction);
    virtual ~PostActEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // Global time of the decay of the primary, as set by the stepping action
    void SetDecayTime(G4double t0) { fDecayTime = t0; }
    G4double GetDecayTime() const  { return fDecayTime; }

private:
    PostActRunAction* fRunAction;
    PostActScintillatorSD* fSD;     // found at the first event
    G4double fDecayTime;
};

#endif
