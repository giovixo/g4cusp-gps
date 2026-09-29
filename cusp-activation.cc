#include <cstdlib>
#include <iostream>
#include <string>

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"

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
#include "SteppingAction.hh"

#include "G4GenericMessenger.hh"


// Usage: cusp-activation [-t nthreads] [macro]
//   no macro     -> interactive (GUI) session
//   -t nthreads  -> number of worker threads (default 1; 0 = all available cores)
static void PrintUsage(const char* prog)
{
    std::cerr << "Usage: " << prog << " [-t nthreads] [macro]" << std::endl
              << "  -t, --threads N   number of worker threads (default 1, 0 = all cores)" << std::endl;
}


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
    
    // Parse the command line
    G4String macroFile;
    G4int nThreads = 1;
    for (G4int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-t" || arg == "--threads")
        {
            char* end = nullptr;
            if (i + 1 < argc) nThreads = std::strtol(argv[++i], &end, 10);
            if (end == nullptr || *end != '\0' || nThreads < 0)
            {
                PrintUsage(argv[0]);
                return 1;
            }
        }
        else if (macroFile.empty() && arg[0] != '-')
        {
            macroFile = arg;
        }
        else
        {
            PrintUsage(argv[0]);
            return 1;
        }
    }

    // Detect interactive mode (if no macro) and define UI session
    G4UIExecutive* ui = 0;
    if (macroFile.empty())
    {
        ui = new G4UIExecutive(argc, argv);
    }

    // Construct the run manager
    // NB: the G4FORCENUMBEROFTHREADS environment variable, if set, overrides -t
#ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    if (nThreads == 0) nThreads = G4Threading::G4GetNumberOfCores();
    runManager->SetNumberOfThreads(nThreads);
#else
    if (nThreads > 1)
        G4cout << "WARNING: Geant4 built without multithreading, -t " << nThreads << " ignored" << G4endl;
    G4RunManager * runManager = new G4RunManager;
#endif

    
    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList());   

    // // User action initialization
    runManager->SetUserInitialization(new UserActionInitialization());
    
    // Initialize G4 kernel (needed before any /gps command in the macros)
    runManager->Initialize();

    
    // Visualization manager construction
    auto visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // Get the pointer to the User Interface manager
    auto uiManager = G4UImanager::GetUIpointer();

    // User commands (/cusp/...)
    auto steppingMessenger = SteppingAction::CreateMessenger();
    
    if (!ui) // Batch mode
    {
        // execute the macro file given on the command line
        G4String command = "/control/execute ";
        uiManager->ApplyCommand(command+macroFile);
        G4cout << "Batch file executed: " << macroFile << G4endl;
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
    delete steppingMessenger;
    delete visManager;
    delete runManager;
    return 0;
}


