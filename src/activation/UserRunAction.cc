#include "UserRunAction.hh"
#include "UserRun.hh"

#include "G4AnalysisManager.hh"

#include <string>


// Constructor: the ntuple is booked once here (not at every run).
// Notice: it must be done the same way in master and workers
UserRunAction::UserRunAction()
{
    fTimer = new G4Timer;

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("csv");
    analysisManager->CreateNtuple("Events", "Events");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleSColumn("Isotope");
    analysisManager->CreateNtupleDColumn("Lifetime");
    analysisManager->CreateNtupleSColumn("Volume");
    analysisManager->FinishNtuple();
}

//UserRunAction::~UserRunAction()
//{   delete fTimer;
//    delete G4AnalysisManager::Instance();
//}

G4Run* UserRunAction::GenerateRun()
{
    return new UserRun;
}


// Defining actions performed at the beginning and/or the end of each run

void UserRunAction::BeginOfRunAction(const G4Run* run)
{
    // Open one output file per run, so that successive /run/beamOn do not overwrite each other
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("scorefile_run" + std::to_string(run->GetRunID()) + ".csv");
    
    if(IsMaster())
    {
        // The run ID is printed at the beginning of each master run
        G4cout << "INFORMATION: Run No. " << run -> GetRunID() << " start." << G4endl;
        fTimer->Start();
//        analysisManager->FillNtupleDColumn(0, run->GetRunID());
//        analysisManager->FillNtupleDColumn(1, run->GetNumberOfEventToBeProcessed());
//        analysisManager->AddNtupleRow();
    }
}


void UserRunAction::EndOfRunAction(const G4Run* run)
{
    // Close and write root file (one for each thread)
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    
    if (IsMaster()){
        fTimer->Stop();
        // The run ID and the number of processed events are printed at the end of each run
        G4int NbOfEvents = run -> GetNumberOfEvent();
        if (NbOfEvents == 0) return;
        G4cout << "INFORMATION: Run No " << run -> GetRunID() << " end." << G4endl;
        G4cout << "INFORMATION: Number of events = " << NbOfEvents << G4endl;
        G4cout << "INFORMATION: Elapsed time = " << *fTimer << G4endl;
    }
}
