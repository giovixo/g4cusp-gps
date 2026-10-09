#include "PostActSteppingAction.hh"
#include "PostActEventAction.hh"
#include "PostActStackingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4HadronicProcessType.hh"


PostActSteppingAction::PostActSteppingAction(PostActEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{;}


PostActSteppingAction::~PostActSteppingAction()
{;}


void PostActSteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();

    if (track->GetParentID() == 0)
    {
        // Primary: at its decay, move the time origin to the decay time. The secondaries
        // of the step are stacked after the stepping action, so the change is effective.
        if (track->GetTrackStatus() != fStopAndKill) return;
        const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (process == nullptr || process->GetProcessType() != fDecay ||
            process->GetProcessSubType() != fRadioactiveDecay) return;

        const G4double t0 = step->GetPostStepPoint()->GetGlobalTime();
        fEventAction->SetDecayTime(t0);
        for (const G4Track* secondary : *step->GetSecondaryInCurrentStep())
        {
            auto s = const_cast<G4Track*>(secondary);
            s->SetGlobalTime(s->GetGlobalTime() - t0);
        }
        return;
    }

    // Long-lived daughter that has just stopped: kill it, keeping the step's energy deposit
    if (track->GetTrackStatus() == fStopButAlive &&
        PostActStackingAction::IsLongLived(track->GetDefinition()))
        track->SetTrackStatus(fStopAndKill);
}
