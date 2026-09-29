#include "UserRun.hh"

UserRun::UserRun()
{
}

UserRun::~UserRun() {;}


void UserRun::RecordEvent(const G4Event* event)
{
    // This method specifies the actions that must be performed at the
    // end of each event (e.g. retrieve information, score, clean up
    // things, etc.)
    
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
