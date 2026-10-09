#include "PostActDetectorConstruction.hh"
#include "PostActScintillatorSD.hh"

#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"

#include <cstdlib>
#include <unordered_map>


PostActDetectorConstruction::PostActDetectorConstruction()
: DetectorConstruction()
{;}


PostActDetectorConstruction::~PostActDetectorConstruction()
{;}


void PostActDetectorConstruction::ConstructSDandField()
{
    // ScintID of each scintillator, from the physical volume names PV-Scatterer_NNN and PV-Absorber_NNN
    std::unordered_map<const G4VPhysicalVolume*, G4int> idMap;
    std::unordered_map<const G4LogicalVolume*, G4int> logicalIds;   // one logical volume per scintillator
    G4int nScatterers = 0;
    G4int nAbsorbers = 0;

    for (G4VPhysicalVolume* pv : *G4PhysicalVolumeStore::GetInstance())
    {
        const G4String& name = pv->GetName();
        G4int id = -1;
        if (name.rfind("PV-Scatterer_", 0) == 0)
        {
            id = std::atoi(name.c_str() + 13) - 1;
            if (id >= 0 && id < 64) ++nScatterers; else id = -1;
        }
        else if (name.rfind("PV-Absorber_", 0) == 0)
        {
            id = std::atoi(name.c_str() + 12) - 1;
            if (id >= 0 && id < 32) { id += 64; ++nAbsorbers; } else id = -1;
        }
        if (id < 0) continue;
        idMap[pv] = id;
        logicalIds[pv->GetLogicalVolume()] = id;
    }

    if (nScatterers != 64 || nAbsorbers != 32 || logicalIds.size() != 96)
    {
        G4ExceptionDescription msg;
        msg << "Expected 64 scatterers and 32 absorbers (one logical volume each), found "
            << nScatterers << ", " << nAbsorbers << " and " << logicalIds.size() << " logical volumes";
        G4Exception("PostActDetectorConstruction::ConstructSDandField", "PostAct010", FatalException, msg);
    }

    // One sensitive detector per thread, attached to the 96 logical volumes
    auto sd = new PostActScintillatorSD("PostActScintillatorSD", idMap);
    G4SDManager::GetSDMpointer()->AddNewDetector(sd);
    for (const auto& item : logicalIds)
        SetSensitiveDetector(const_cast<G4LogicalVolume*>(item.first), sd);
}
