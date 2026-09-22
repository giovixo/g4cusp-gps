#include <iostream>
#include <map>

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "DetectorConstruction.hh"
#include "UserActionInitialization.hh"

#include "PhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "UserRunAction.hh"
#include "UserEventAction.hh"

#include "G4GenericMessenger.hh"


#include "GlobalRsmSource.hh"

int main(int argc, char **argv)
{
    // Detect the C++ standard version (C++17 is the GEANT4 recomandation)
    std::cout << "C++ Standard Version: ";
    #if __cplusplus == 201703L
       std::cout << "C++17" << std::endl;
    #elif __cplusplus == 202002L
       std::cout << "C++20" << std::endl;
    #elif __cplusplus == 201402L
       std::cout << "C++14" << std::endl;
    #elif __cplusplus == 201103L
       std::cout << "C++11" << std::endl;
    #else
       std::cout << "Pre-C++11 or unknown standard" << std::endl;
    #endif

//    G4Random::setTheEngine(new CLHEP::RanecuEngine);
//    G4Random::setTheSeed(time(0));
    
    // Detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = 0;
    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    // Construct the run manager
#ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(1);
#else
    G4RunManager * runManager = new G4RunManager;
#endif

    // Initialize the global RsmSource
    gRsmSource = new RsmSource();
    gRsmSource->ReadCSV("rsm-events.csv");
    #ifdef DEBUG
       gRsmSource->PrintCSV(); // Print the CSV data for debugging
    #endif
    
    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList());   

    // // User action initialization
    runManager->SetUserInitialization(new UserActionInitialization());
    
    // Initialize G4 kernel
    runManager->Initialize();

    
    // Visualization manager construction
    auto visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // Get the pointer to the User Interface manager
    auto uiManager = G4UImanager::GetUIpointer();
    
    uiManager->ApplyCommand("/run/initialize");
    
    if (!ui) // Batch mode
    {
        // execute an argument macro file if exist
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        uiManager->ApplyCommand(command+fileName);
        G4cout << "Batch file executed: " << fileName << G4endl;
    }
    else // GUI or interactive mode
    {
        uiManager->ApplyCommand("/control/execute init_vis.mac");
        
        if (ui->IsGUI()) {
            uiManager->ApplyCommand("/control/execute visGUI.mac");
            uiManager->ApplyCommand("/control/execute visGUI2.mac");
        }
        // start interactive session
        ui->SessionStart();
        delete ui;
    }

    
    // Job termination
    delete visManager;
    delete runManager;
    return 0;
}


