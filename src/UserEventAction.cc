#include "UserRunAction.hh"
#include "UserEventAction.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>
#include <vector>

// Defining actions performed at the beginning and/or the end of each event

// Constructor
UserEventAction::UserEventAction()
{;}

// Destructor
UserEventAction::~UserEventAction()
{;}



void UserEventAction::BeginOfEventAction(const G4Event* event) {
    // This method specifies the actions that must be performed at the
    // beginning of each event.
    G4int nEvent = event -> GetEventID()+1;
    if(!(nEvent % 10000))
    {
        G4cout << "INFORMATION: event number " << nEvent << G4endl;
    }
    
}


void UserEventAction::EndOfEventAction(const G4Event* event) {
    // This method specifies the actions that must be performed at the
    // end of each event.
    //G4cout << "---> UserEventAction: End of event: " << event -> GetEventID() << G4endl;
    
    // The number of events is retrieved from the G4Event object
    // Remember: the first event has number 0, therefore we add 1
    G4int nmbEvents = (event -> GetEventID()) + 1;
    
    // The average CPU time is printed if a run has at least 10000 events
    if(!(nmbEvents % 100000))
    {
        G4cout << "INFORMATION: " << nmbEvents << " events processed." << G4endl;
    }
    
}
