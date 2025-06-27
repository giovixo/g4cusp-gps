#ifndef USERRUNACTION_HH
#define USERRUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "globals.hh"

#include "G4UImanager.hh"
#include "G4GenericMessenger.hh"

#include "G4AnalysisManager.hh"

#include "G4UIcommand.hh"

// Defining actions performed at the beginning and/or the end of each run

// G4UserRunAction is the base class for defining user actions performed at
// the beginning and/or the end of each run.


class UserRunAction : public G4UserRunAction
{
public:
    UserRunAction() {
    //fOutFileName = "scorefile.root";  
    // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    fTimer = new G4Timer;
    }     // Constructor
    ~UserRunAction() {delete fTimer;}  // Destructor

    void SetFileName(G4String file_name) {fOutFileName = file_name;}
    G4String GetFileName() {
        return fOutFileName;
    }
    
    // G4UserRunAction has two methods, BeginOfRunAction and EndOfRunAction,
    // which can be overloaded by the user to define specific actions performed
    // the beginning and the end of each run.
    // Information about the run can be retrieved from the G4Run objects passed
    // to these functions.
    G4Run* GenerateRun();
    void BeginOfRunAction(const G4Run* run);
    void   EndOfRunAction(const G4Run* run);
       
private:
    G4Timer* fTimer;
    G4String fOutFileName;    
};

#endif // USERRUNACTION_HH
