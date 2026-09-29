#include "UserRun.hh"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcommand.hh"
#include "tls.hh"

#include "ConfigFile.hh"

//#include "g4root.hh"
#include "G4AnalysisManager.hh"

UserRun::UserRun()
{
}

UserRun::~UserRun() {;}


void UserRun::RecordEvent(const G4Event* event)
{
    // This method specifies the actions that must be performed at the
    // end of each event (e.g. retrieve information, score, clean up
    // things, etc.)
    
    G4double eventID = event -> GetEventID();
    // Record the event
    //G4cout << "---> (Record Event) End of event: " << event -> GetEventID() << G4endl;
    G4Run::RecordEvent(event);
}

void UserRun::Merge(const G4Run* run)
{
    G4cout << "---> Merging files... " << G4endl;
    G4Run::Merge(run);
    G4cout << "---> Merged files. " << G4endl;

}
