#include "PostActRunAction.hh"
#include "PostActConfig.hh"

#include "G4HadronicParameters.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"

#include <filesystem>
#include <iomanip>


PostActRunAction::PostActRunAction()
: G4UserRunAction(), fTimer(new G4Timer)
{;}


PostActRunAction::~PostActRunAction()
{
    if (fFile.is_open()) fFile.close();
    delete fTimer;
}


// Opens (appending) a file; writes the header if the file is new or empty
void PostActRunAction::OpenFile(const G4String& fileName, const G4String& header)
{
    if (fFile.is_open()) fFile.close();
    std::error_code ec;
    const std::string path(fileName);
    const bool empty = !std::filesystem::exists(path, ec) || std::filesystem::file_size(path, ec) == 0;
    fFile.open(fileName, std::ios::out | std::ios::app);
    if (!fFile.is_open())
    {
        G4ExceptionDescription msg;
        msg << "Cannot open " << fileName;
        G4Exception("PostActRunAction::OpenFile", "PostAct012", FatalException, msg);
    }
    if (empty) fFile << header << "\n";
    fFileName = fileName;
}


void PostActRunAction::BeginOfRunAction(const G4Run*)
{
    if (IsMaster())
    {
        // Decays after the threshold are ignored by Geant4 (default 1 year!). It is set to
        // 1e60 y in main, before the physics is constructed; this only checks it.
        const G4double threshold = G4HadronicParameters::Instance()->GetTimeThresholdForRadioactiveDecay();
        if (threshold < 1.e50*year)
            G4cout << "WARNING: decays later than " << G4BestUnit(threshold, "Time")
                   << " are ignored by the radioactive decay" << G4endl;
        fTimer->Start();
        return;
    }

    // Worker: one file per thread, opened at the first run (or if the prefix changed)
    const G4int threadID = G4Threading::IsMultithreadedApplication() ? G4Threading::G4GetThreadId() : 0;
    const G4String fileName = PostActConfig::OutputPrefix() + "_t" + std::to_string(threadID) + ".csv";
    if (fileName != fFileName || !fFile.is_open())
        OpenFile(fileName, "RunID,EventID,ScintID,Edep_keV,t_ns");
}


void PostActRunAction::EndOfRunAction(const G4Run* run)
{
    if (!IsMaster())
    {
        if (fFile.is_open()) fFile.flush();
        return;
    }

    fTimer->Stop();
    const G4int nDecays = run->GetNumberOfEvent();
    const G4double elapsed = fTimer->GetRealElapsed();

    OpenFile(PostActConfig::OutputPrefix() + "_runs.csv", "RunID,Volume,Isotope,NDecays");
    fFile << run->GetRunID() << "," << PostActConfig::VolumeName() << ","
          << PostActConfig::IsotopeName() << "," << nDecays << "\n";
    fFile.close();

    G4cout << "INFORMATION: run " << run->GetRunID() << ": " << nDecays << " decays of "
           << PostActConfig::IsotopeName() << " in " << PostActConfig::VolumeName()
           << ", " << elapsed << " s";
    if (elapsed > 0.) G4cout << " (" << nDecays/elapsed << " decays/s)";
    G4cout << G4endl;
}


void PostActRunAction::WriteRow(G4int runID, G4int eventID, G4int scintID, G4double edep, G4double time)
{
    fFile << runID << "," << eventID << "," << scintID << ","
          << std::setprecision(6) << edep/keV << "," << time/ns << "\n";
}
