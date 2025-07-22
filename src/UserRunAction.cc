#include "UserRunAction.hh"
#include "G4Run.hh"
#include "UserRun.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "ConfigFile.hh"
// Change this header if you want different AnalysisManager output (e.g. XML)
//#include "g4root.hh"
//#include "G4AnalysisManager.hh"

#include <fstream>
#include <vector>
#include <string>


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
    // Create analysis manager
    // Notice: it must be done the same way in master and workers
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);


    // Open the output file name
    G4String fileName = analysisManager->GetFileName();
    if (!fileName.empty()) {
        G4cout << "Output file: " << fileName << G4endl;
        analysisManager->OpenFile(fileName);
    } else {
        G4cout << "Output file: " << "scorefile.root" << G4endl;
        analysisManager->OpenFile("scorefile.root");
    }

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstNtupleId(1);
    // Ntuple merging (only for Geant v. 4.10.03 or higher)
    analysisManager->SetNtupleMerging(true);
    
    analysisManager->CreateNtuple("Events", "Events");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleDColumn("Timestamp");
    analysisManager->CreateNtupleDColumn("En_dep");
    analysisManager->CreateNtupleIColumn("Scint_ID");
    analysisManager->CreateNtupleDColumn("X_Primary");
    analysisManager->CreateNtupleDColumn("Y_Primary");
    analysisManager->CreateNtupleDColumn("Z_Primary");
    analysisManager->CreateNtupleDColumn("Theta_Primary");
    analysisManager->CreateNtupleDColumn("Phi_Primary");
    analysisManager->CreateNtupleDColumn("En_Primary");
    analysisManager->CreateNtupleDColumn("Event_time");
    analysisManager->CreateNtupleDColumn("X_Detected");
    analysisManager->CreateNtupleDColumn("Y_Detected");
    analysisManager->CreateNtupleDColumn("Z_Detected");
    analysisManager->CreateNtupleDColumn("X_Pol");
    analysisManager->CreateNtupleDColumn("Y_Pol");
    analysisManager->CreateNtupleDColumn("Z_Pol");
    analysisManager->FinishNtuple();
    
    
    if(IsMaster())
    {
        // The run ID is printed at the beginning of each master run
        G4cout << "INFORMATION: Run No. " << run -> GetRunID() << " start." << G4endl;
        fTimer->Start();
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
        G4cout << "INFORMATION: Elapsed time: " << *fTimer << G4endl;
    }
}
