#ifndef PostActRunAction_h
#define PostActRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <fstream>
#include <set>

class G4Timer;

// Workers: write the rows of the events to <prefix>_t<threadID>.csv
//          (RunID,EventID,ScintID,Edep_keV,t_ns), opened at the first run.
// Master:  appends one row per run to <prefix>_runs.csv (RunID,Volume,Isotope,NDecays).
// Existing files are appended to (the header is written only in a new or empty file).

class PostActRunAction : public G4UserRunAction
{
public:
    PostActRunAction();
    virtual ~PostActRunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

    // Worker only
    void WriteRow(G4int runID, G4int eventID, G4int scintID, G4double edep, G4double time);

private:
    void OpenFile(const G4String& fileName, const G4String& header);

    std::ofstream fFile;
    G4String fFileName;
    G4Timer* fTimer;
};

#endif
