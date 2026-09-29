#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

SteppingAction::SteppingAction()
{ }


SteppingAction::~SteppingAction()
{ }


void SteppingAction::UserSteppingAction(const G4Step* step)
{

      // Which volume?

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();

      G4String volumeName = preStepPoint->GetTouchableHandle()->GetVolume()->GetName();

    const G4VProcess* currentProcess0 = preStepPoint->GetProcessDefinedStep();
    const G4VProcess* currentProcess1 = postStepPoint->GetProcessDefinedStep();
    const G4VProcess* creatorProcess = step->GetTrack()->GetCreatorProcess();

    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
    G4double globaltime = step->GetTrack()-> GetGlobalTime();
    G4int eventID;

    if (!(particle->GetPDGStable()) && !(particleName == "neutron"))
    {
        if (particle->GetPDGLifeTime() > 1e-1*s && particle->GetPDGLifeTime() < 1e18*s) {
            
            eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
            G4cout << "*** RADIOISOTOPE " << eventID;
            G4cout << " " << particleName;
            G4cout << " " << particle->GetAtomicNumber();
            G4cout << " " << particle->GetAtomicMass();
            G4cout << " " << particle->GetPDGEncoding();
            G4cout << " " << particle->GetPDGLifeTime()/s;
            G4cout << " " << volumeName << G4endl;
            
            
            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            analysisManager -> FillNtupleIColumn(0, eventID);
            analysisManager -> FillNtupleSColumn(1, particleName);
            analysisManager -> FillNtupleDColumn(2, particle->GetPDGLifeTime()/s);
            analysisManager -> FillNtupleSColumn(3, volumeName);
            analysisManager -> AddNtupleRow();


//            if (currentProcess0 != 0) G4cout << "*** DEBUG - POSSIBLE RADIOACTIVE DECAY pre "  << currentProcess0->GetProcessName() << G4endl;
//            if (currentProcess1 != 0) G4cout << "*** DEBUG - POSSIBLE RADIOACTIVE DECAY post " << currentProcess1->GetProcessName() << G4endl;
//            if (creatorProcess != 0) G4cout << "*** DEBUG - POSSIBLE RADIOACTIVE DECAY creator " << creatorProcess->GetProcessName() << G4endl;
    //        particle->DumpTable();
//            G4cout << G4endl;

            step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
        }
        
    }

  
}

