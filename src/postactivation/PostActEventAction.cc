#include "PostActEventAction.hh"
#include "PostActRunAction.hh"
#include "PostActScintillatorSD.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"


PostActEventAction::PostActEventAction(PostActRunAction* runAction)
: G4UserEventAction(), fRunAction(runAction), fSD(nullptr), fDecayTime(0.)
{;}


PostActEventAction::~PostActEventAction()
{;}


void PostActEventAction::BeginOfEventAction(const G4Event*)
{
    fDecayTime = 0.;
    if (fSD == nullptr)
    {
        // The sensitive detector is thread-local: G4SDManager returns the one of this thread
        fSD = dynamic_cast<PostActScintillatorSD*>(
                G4SDManager::GetSDMpointer()->FindSensitiveDetector("PostActScintillatorSD", false));
        if (fSD == nullptr)
            G4Exception("PostActEventAction::BeginOfEventAction", "PostAct020", FatalException,
                        "Sensitive detector PostActScintillatorSD not found");
    }
}


void PostActEventAction::EndOfEventAction(const G4Event* event)
{
    const G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    for (G4int id = 0; id < PostActScintillatorSD::kNScintillators; ++id)
    {
        const G4double edep = fSD->GetEdep(id);
        if (edep > 0.)
            fRunAction->WriteRow(runID, event->GetEventID(), id, edep, fSD->GetTime(id));
    }
    fSD->Reset();
}
