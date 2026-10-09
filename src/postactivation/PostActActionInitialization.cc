#include "PostActActionInitialization.hh"
#include "PostActPrimaryGeneratorAction.hh"


PostActActionInitialization::PostActActionInitialization()
: G4VUserActionInitialization()
{;}


PostActActionInitialization::~PostActActionInitialization()
{;}


void PostActActionInitialization::BuildForMaster() const
{
}

void PostActActionInitialization::Build() const
{
    SetUserAction(new PostActPrimaryGeneratorAction);

    // TODO: stepping action (kill daughters with T1/2 >= 1 us), sensitive detectors and output
}
