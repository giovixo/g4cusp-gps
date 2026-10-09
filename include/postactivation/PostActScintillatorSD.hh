#ifndef PostActScintillatorSD_h
#define PostActScintillatorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Types.hh"

#include <array>
#include <unordered_map>

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;
class G4VPhysicalVolume;

// Sensitive detector attached to the 96 scintillators (64 scatterers + 32 absorbers).
// It keeps, for each scintillator, the energy deposited in the event and the time of
// the first deposit. Only deposits within PostActConfig::TimeWindow() after the decay
// are counted (the time origin is the decay, see PostActSteppingAction).
// The totals are a fixed array, reset after each event: no hits collection.
// One instance per thread; the event action finds it by name through G4SDManager.

class PostActScintillatorSD : public G4VSensitiveDetector
{
public:
    static constexpr G4int kNScintillators = 96;   // IDs 0-63: scatterers, 64-95: absorbers

    PostActScintillatorSD(const G4String& name,
                          const std::unordered_map<const G4VPhysicalVolume*, G4int>& idMap);
    virtual ~PostActScintillatorSD();

    virtual void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*);

    // Zeroes the totals (called by the event action at the end of each event)
    void Reset();

    G4double GetEdep(G4int id) const { return fEdep[id]; }
    G4double GetTime(G4int id) const { return fTime[id]; }

private:
    std::unordered_map<const G4VPhysicalVolume*, G4int> fIdMap;
    std::array<G4double, kNScintillators> fEdep;
    std::array<G4double, kNScintillators> fTime;
    G4double fTimeWindow;
};

#endif
