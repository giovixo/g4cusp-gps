#include "PostActScintillatorSD.hh"
#include "PostActConfig.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <limits>


PostActScintillatorSD::PostActScintillatorSD(const G4String& name,
        const std::unordered_map<const G4VPhysicalVolume*, G4int>& idMap)
: G4VSensitiveDetector(name), fIdMap(idMap), fTimeWindow(PostActConfig::TimeWindow())
{
    Reset();
}


PostActScintillatorSD::~PostActScintillatorSD()
{;}


void PostActScintillatorSD::Reset()
{
    fEdep.fill(0.);
    fTime.fill(std::numeric_limits<G4double>::max());
}


void PostActScintillatorSD::Initialize(G4HCofThisEvent*)
{
    Reset();
    fTimeWindow = PostActConfig::TimeWindow();   // may change between runs
}


G4bool PostActScintillatorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return false;

    // Time of the deposit: start of the step (the time origin is the decay)
    const G4double time = step->GetPreStepPoint()->GetGlobalTime();
    if (time > fTimeWindow) return false;

    auto it = fIdMap.find(step->GetPreStepPoint()->GetPhysicalVolume());
    if (it == fIdMap.end()) return false;

    fEdep[it->second] += edep;
    fTime[it->second] = std::min(fTime[it->second], time);
    return true;
}
