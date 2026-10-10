#include "PostActRunAction.hh"
#include "PostActConfig.hh"

#include "G4HadronicParameters.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"

#include <filesystem>
#include <iomanip>
#include <regex>


namespace { G4Mutex openMutex = G4MUTEX_INITIALIZER; }

std::set<G4String> PostActRunAction::fOpenedFiles;


PostActRunAction::PostActRunAction()
: G4UserRunAction(), fTimer(new G4Timer)
{;}


PostActRunAction::~PostActRunAction()
{
    if (fFile.is_open()) fFile.close();
    delete fTimer;
}


// Deletes the output files of a previous job with this prefix (also thread files that this
// job would not overwrite, e.g. if it runs with fewer threads)
void PostActRunAction::RemoveOldFiles(const G4String& prefix)
{
    namespace fs = std::filesystem;
    const fs::path base(prefix.c_str());
    const fs::path dir = base.has_parent_path() ? base.parent_path() : fs::path(".");
    const std::regex pattern(std::regex_replace(base.filename().string(),
                                                std::regex(R"([.^$|()\[\]{}*+?\\])"), R"(\$&)")
                             + R"(_(t\d+|runs)\.csv)");
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(dir, ec))
        if (std::regex_match(entry.path().filename().string(), pattern))
            fs::remove(entry.path(), ec);
}


// Opens a file: truncated (with the header) the first time this process opens it, appended afterwards
void PostActRunAction::OpenFile(const G4String& fileName, const G4String& header)
{
    if (fFile.is_open()) fFile.close();
    G4bool first;
    {
        G4AutoLock lock(&openMutex);
        first = fOpenedFiles.insert(fileName).second;
    }
    fFile.open(fileName, std::ios::out | (first ? std::ios::trunc : std::ios::app));
    if (!fFile.is_open())
    {
        G4ExceptionDescription msg;
        msg << "Cannot open " << fileName;
        G4Exception("PostActRunAction::OpenFile", "PostAct021", FatalException, msg);
    }
    if (first) fFile << header << "\n";
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
        // First run of this process with this prefix: remove the files of an earlier job.
        // The master runs this before the workers start the run.
        if (fOpenedFiles.count(PostActConfig::OutputPrefix() + "_runs.csv") == 0)
            RemoveOldFiles(PostActConfig::OutputPrefix());
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

    OpenFile(PostActConfig::OutputPrefix() + "_runs.csv", "RunID,Volume,Isotope,NDecays,G4Ion");
    fFile << run->GetRunID() << "," << PostActConfig::VolumeName() << ","
          << PostActConfig::IsotopeName() << "," << nDecays << "," << PostActConfig::G4IonName() << "\n";
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
