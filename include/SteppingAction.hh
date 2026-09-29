#ifndef USERSTEPPINGACTION_HH
#define USERSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4GenericMessenger;

// Records every radioactive nuclide produced in the geometry (ntuple "Events")
// and kills it together with its secondaries, so that it is not decayed.

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
   ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

    // Creates the /cusp/stepping/ commands. Call it once, on the master thread (see main):
    // SteppingAction objects live only on the worker threads.
    static G4GenericMessenger* CreateMessenger();

  private:
    // Lifetime window of the recorded nuclides
    static constexpr G4double fMinLifetime = 0.1*s;
    static constexpr G4double fMaxLifetime = 1e18*s;

    // 0 = silent, 1 = one "*** RADIOISOTOPE" line per recorded nuclide.
    // Shared by all threads; set by /cusp/stepping/verbose before /run/beamOn.
    static G4int fVerbose;
};

#endif
