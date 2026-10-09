#include "PostActActionInitialization.hh"
#include "PostActPrimaryGeneratorAction.hh"
#include "PostActRunAction.hh"
#include "PostActEventAction.hh"
#include "PostActSteppingAction.hh"
#include "PostActStackingAction.hh"


PostActActionInitialization::PostActActionInitialization()
: G4VUserActionInitialization()
{;}


PostActActionInitialization::~PostActActionInitialization()
{;}


void PostActActionInitialization::BuildForMaster() const
{
    SetUserAction(new PostActRunAction);
}

void PostActActionInitialization::Build() const
{
    SetUserAction(new PostActPrimaryGeneratorAction);

    auto runAction = new PostActRunAction;
    SetUserAction(runAction);
    auto eventAction = new PostActEventAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(new PostActSteppingAction(eventAction));
    SetUserAction(new PostActStackingAction);
}
