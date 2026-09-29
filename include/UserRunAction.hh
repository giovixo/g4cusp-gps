#ifndef USERRUNACTION_HH
#define USERRUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "globals.hh"

// Defining actions performed at the beginning and/or the end of each run

// G4UserRunAction is the base class for defining user actions performed at
// the beginning and/or the end of each run.

class UserRunAction : public G4UserRunAction
{
public:
    UserRunAction() {fTimer = new G4Timer;}        // Constructor
    ~UserRunAction() {delete fTimer;}  // Destructor
    
    // G4UserRunAction has two methods, BeginOfRunAction and EndOfRunAction,
    // which can be overloaded by the user to define specific actions performed
    // the beginning and the end of each run.
    // Information about the run can be retrieved from the G4Run objects passed
    // to these functions.
    G4Run* GenerateRun();
    void BeginOfRunAction(const G4Run* run);
    void   EndOfRunAction(const G4Run* run);
    
    
private:
    G4Timer* fTimer;
    
};

#endif // USERRUNACTION_HH
