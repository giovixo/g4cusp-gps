#include "UserActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "UserRunAction.hh"
#include "UserEventAction.hh"


UserActionInitialization::UserActionInitialization()
: G4VUserActionInitialization()
{;}


UserActionInitialization::~UserActionInitialization()
{;}


void UserActionInitialization::BuildForMaster() const
{
    SetUserAction(new UserRunAction);
}

void UserActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new UserRunAction);
    SetUserAction(new UserEventAction);
}
