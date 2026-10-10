#include "PostActDetectorConstruction.hh"
#include "PostActScintillatorSD.hh"

#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Material.hh"
#include "G4Threading.hh"

#include <array>
#include <regex>
#include <unordered_map>


PostActDetectorConstruction::PostActDetectorConstruction()
: DetectorConstruction()
{;}


PostActDetectorConstruction::~PostActDetectorConstruction()
{;}


void PostActDetectorConstruction::ConstructSDandField()
{
    // ScintID from the physical volume name, matched exactly:
    //   PV-Scatterer_NNN (plastic) -> NNN-1       (0-63)
    //   PV-Absorber_NNN  (GAGG)    -> 64 + NNN-1  (64-95)
    // Other volumes with these prefixes (PV-Scatterer_Filter_Ti_001, PV-Absorber_Collimator_002, ...)
    // do not match. Each ID must be found exactly once, in the expected material.
    static const std::regex pattern(R"(^PV-(Scatterer|Absorber)_(\d{3})$)");
    std::unordered_map<const G4VPhysicalVolume*, G4int> idMap;
    std::unordered_map<const G4LogicalVolume*, G4int> logicalIds;   // one logical volume per scintillator
    std::array<const G4VPhysicalVolume*, PostActScintillatorSD::kNScintillators> byId{};

    G4ExceptionDescription errors;
    for (G4VPhysicalVolume* pv : *G4PhysicalVolumeStore::GetInstance())
    {
        std::smatch m;
        const std::string name = pv->GetName();
        if (!std::regex_match(name, m, pattern)) continue;

        const G4bool scatterer = (m[1] == "Scatterer");
        const G4int n = std::stoi(m[2].str());
        const G4int nMax = scatterer ? 64 : 32;
        const G4String material = pv->GetLogicalVolume()->GetMaterial()->GetName();
        const G4String expected = scatterer ? "G4_PLASTIC_SC_VINYLTOLUENE" : "GAGG";
        if (n < 1 || n > nMax)
        {
            errors << "  " << name << ": number out of range 1-" << nMax << "\n";
            continue;
        }
        if (material != expected)
            errors << "  " << name << ": material " << material << ", expected " << expected << "\n";

        const G4int id = (scatterer ? 0 : 64) + n - 1;
        if (byId[id] != nullptr)
            errors << "  " << name << ": ScintID " << id << " already used by " << byId[id]->GetName() << "\n";
        byId[id] = pv;
        idMap[pv] = id;
        logicalIds[pv->GetLogicalVolume()] = id;
    }
    for (G4int id = 0; id < PostActScintillatorSD::kNScintillators; ++id)
        if (byId[id] == nullptr) errors << "  ScintID " << id << ": no volume\n";
    if (logicalIds.size() != idMap.size())
        errors << "  some scintillators share a logical volume\n";

    if (!errors.str().empty())
    {
        G4ExceptionDescription msg;
        msg << "Scintillator mapping (64 PV-Scatterer_NNN + 32 PV-Absorber_NNN) failed:\n" << errors.str();
        G4Exception("PostActDetectorConstruction::ConstructSDandField", "PostAct010", FatalException, msg);
    }
    if (G4Threading::G4GetThreadId() <= 0)   // once: first worker, or sequential mode
        G4cout << "PostAct: ScintID 0-63 = " << byId[0]->GetName() << " ... " << byId[63]->GetName()
               << ", 64-95 = " << byId[64]->GetName() << " ... " << byId[95]->GetName() << G4endl;

    // One sensitive detector per thread, attached to the 96 logical volumes
    auto sd = new PostActScintillatorSD("PostActScintillatorSD", idMap);
    G4SDManager::GetSDMpointer()->AddNewDetector(sd);
    for (const auto& item : logicalIds)
        SetSensitiveDetector(const_cast<G4LogicalVolume*>(item.first), sd);
}
