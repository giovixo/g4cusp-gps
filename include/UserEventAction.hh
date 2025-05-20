#ifndef USEREVENTACTION_HH
#define USEREVENTACTION_HH

#include "G4Timer.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"


// Defining actions performed at the beginning and/or the end of each event

// G4UserEventAction is the base class for defining user actions performed at
// the beginning and/or end of each event.
class UserRunAction;

class UserEventAction : public G4UserEventAction
{
public:
    UserEventAction();	// Constructor
    virtual ~UserEventAction();		    // Destructor
    
    // G4UserEventAction has two methods, BeginOfEventAction and
    // EndOfEventAction, which can be overloaded by the user to define specific
    // actions performed at the beginning and at the end of each event.
    // Information about the event can be retrieved from the G4Event objects
    // passed to these functions.
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    
};

#endif // USEREVENTACTION_HH
