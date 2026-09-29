#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4Neutron.hh"
#include "G4GenericMessenger.hh"
#include "G4AnalysisManager.hh"

G4int SteppingAction::fVerbose = 1;

SteppingAction::SteppingAction()
{ }


SteppingAction::~SteppingAction()
{ }


G4GenericMessenger* SteppingAction::CreateMessenger()
{
    auto messenger = new G4GenericMessenger(nullptr, "/cusp/stepping/", "Radioactive nuclide recording");
    messenger->DeclareProperty("verbose", fVerbose,
                               "0 = silent, 1 = print one line per recorded nuclide")
        .SetParameterName("level", false)
        .SetRange("level>=0")
        .SetDefaultValue("1")
        .SetToBeBroadcasted(false);   // the value is shared: set it on the master only
    return messenger;
}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    const G4ParticleDefinition* particle = track->GetDefinition();

    // Radioactive nuclides: unstable (but not neutrons) with a lifetime in the window
    if (particle->GetPDGStable() || particle == G4Neutron::Definition()) return;

    const G4double lifetime = particle->GetPDGLifeTime();
    if (lifetime <= fMinLifetime || lifetime >= fMaxLifetime) return;

    // Volume where the nuclide is (its first step starts where it was produced)
    const G4String& volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
    const G4String& particleName = particle->GetParticleName();
    const G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    if (fVerbose > 0)
    {
        G4cout << "*** RADIOISOTOPE " << eventID
               << " " << particleName
               << " " << particle->GetAtomicNumber()
               << " " << particle->GetAtomicMass()
               << " " << particle->GetPDGEncoding()
               << " " << lifetime/s
               << " " << volumeName << G4endl;
    }

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(0, eventID);
    analysisManager->FillNtupleSColumn(1, particleName);
    analysisManager->FillNtupleDColumn(2, lifetime/s);
    analysisManager->FillNtupleSColumn(3, volumeName);
    analysisManager->AddNtupleRow();

    // Do not decay it: kill the nuclide and anything produced in this step
    track->SetTrackStatus(fKillTrackAndSecondaries);
}
