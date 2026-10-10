#include "PostActPrimaryGeneratorAction.hh"
#include "PostActConfig.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4ParticleGun.hh"
#include "G4Point3D.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"

namespace
{
    const G4int kMaxTries = 1000000;            // maximum tries to find a point inside the solid
    const G4int kTestPoints = 10000;            // points used to measure the acceptance of a new volume
}


PostActPrimaryGeneratorAction::PostActPrimaryGeneratorAction()
    : fIon(nullptr), fSolid(nullptr)
{
    // The gun only carries the particle at rest; the position is set at every event
    fParticleGun = new G4ParticleGun(1);
    fParticleGun->SetParticleEnergy(0.);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}


PostActPrimaryGeneratorAction::~PostActPrimaryGeneratorAction()
{
    delete fParticleGun;
}


void PostActPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    UpdateIon();
    UpdateVolume();

    // Rejection sampling in the bounding box; no fallback if the solid is never hit
    G4ThreeVector local;
    if (!Sample(local))
    {
        G4ExceptionDescription msg;
        msg << "No point inside \"" << fCurrentVolume << "\" found in " << kMaxTries << " tries";
        G4Exception("PostActPrimaryGeneratorAction::GeneratePrimaries", "PostAct014",
                    FatalException, msg);
    }

    fParticleGun->SetParticlePosition(fTransform * G4Point3D(local));
    fParticleGun->GeneratePrimaryVertex(anEvent);
}


// Uniform point inside the solid (local frame); false if kMaxTries are not enough
G4bool PostActPrimaryGeneratorAction::Sample(G4ThreeVector& local) const
{
    for (G4int i = 0; i < kMaxTries; ++i)
    {
        // Each axis has its own limits
        local.set(fLo.x() + (fHi.x() - fLo.x())*G4UniformRand(),
                  fLo.y() + (fHi.y() - fLo.y())*G4UniformRand(),
                  fLo.z() + (fHi.z() - fLo.z())*G4UniformRand());
        if (fSolid->Inside(local) == kInside) return true;
    }
    return false;
}


void PostActPrimaryGeneratorAction::UpdateIon()
{
    const G4String& name = PostActConfig::IsotopeName();
    if (name.empty())
    {
        G4Exception("PostActPrimaryGeneratorAction::UpdateIon", "PostAct011", FatalException,
                    "Isotope not set: use /postact/isotope <name> before /run/beamOn");
        return;
    }
    if (fIon != nullptr && name == fCurrentIsotope) return;

    const PostActIon& req = PostActConfig::Ion();
    const G4Ions::G4FloatLevelBase flb = (req.flb == '\0') ? G4Ions::G4FloatLevelBase::no_Float
                                                           : G4Ions::FloatLevelBase(req.flb);
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(req.Z, req.A, req.E, flb);

    // Same ion as created on the master by /postact/isotope (see PostActConfig::SetIsotope)
    if (ion == nullptr || ion->GetParticleName() != PostActConfig::G4IonName())
    {
        G4ExceptionDescription msg;
        msg << "Isotope " << name << ": the worker got "
            << (ion == nullptr ? G4String("no ion") : ion->GetParticleName())
            << " instead of " << PostActConfig::G4IonName();
        G4Exception("PostActPrimaryGeneratorAction::UpdateIon", "PostAct013", FatalException, msg);
        return;
    }

    fIon = ion;
    fCurrentIsotope = name;
    fParticleGun->SetParticleDefinition(fIon);
    fParticleGun->SetParticleCharge(0.);    // neutral at rest
    fParticleGun->SetParticleEnergy(0.);

    G4cout << "PostAct generator: " << fIon->GetParticleName() << ", ";
    if (fIon->GetPDGStable()) G4cout << "stable";
    else G4cout << "unstable, lifetime " << G4BestUnit(fIon->GetPDGLifeTime(), "Time");
    G4cout << G4endl;
}


// Placement chain from mother down to target; transform = local(target) -> frame of mother
G4bool PostActPrimaryGeneratorAction::FindPath(const G4VPhysicalVolume* mother,
                                               const G4VPhysicalVolume* target,
                                               G4Transform3D& transform)
{
    const G4LogicalVolume* lv = mother->GetLogicalVolume();
    for (std::size_t i = 0; i < lv->GetNoDaughters(); ++i)
    {
        const G4VPhysicalVolume* d = lv->GetDaughter(i);
        // Object rotation = inverse of the frame rotation stored in the PV
        G4Transform3D t(d->GetObjectRotationValue(), d->GetTranslation());
        if (d == target)
        {
            transform = t;
            return true;
        }
        G4Transform3D sub;
        if (FindPath(d, target, sub))
        {
            transform = t * sub;
            return true;
        }
    }
    return false;
}


void PostActPrimaryGeneratorAction::UpdateVolume()
{
    const G4String& name = PostActConfig::VolumeName();
    if (name.empty())
    {
        G4Exception("PostActPrimaryGeneratorAction::UpdateVolume", "PostAct015", FatalException,
                    "Volume not set: use /postact/volume <PV name> before /run/beamOn");
        return;
    }
    if (fSolid != nullptr && name == fCurrentVolume) return;

    G4ExceptionDescription msg;
    msg << "Volume " << name << ": ";
    const G4VPhysicalVolume* pv = G4PhysicalVolumeStore::GetInstance()->GetVolume(name, false);
    if (pv == nullptr)
    {
        msg << "not found";
        G4Exception("PostActPrimaryGeneratorAction::UpdateVolume", "PostAct016", FatalException, msg);
        return;
    }
    if (pv->IsReplicated() || pv->IsParameterised())
    {
        msg << "replicated or parameterised volumes are not supported";
        G4Exception("PostActPrimaryGeneratorAction::UpdateVolume", "PostAct017", FatalException, msg);
        return;
    }

    // Local -> global transform, composed along the path from the world volume
    const G4VPhysicalVolume* world =
        G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
    G4Transform3D transform;
    if (pv != world && !FindPath(world, pv, transform))
    {
        msg << "not found in the geometry tree of the world volume";
        G4Exception("PostActPrimaryGeneratorAction::UpdateVolume", "PostAct018", FatalException, msg);
        return;
    }

    fSolid = pv->GetLogicalVolume()->GetSolid();
    fTransform = transform;
    fCurrentVolume = name;
    fSolid->BoundingLimits(fLo, fHi);

    // Acceptance of the rejection sampling (information only)
    G4int accepted = 0;
    for (G4int i = 0; i < kTestPoints; ++i)
    {
        const G4ThreeVector p(fLo.x() + (fHi.x() - fLo.x())*G4UniformRand(),
                              fLo.y() + (fHi.y() - fLo.y())*G4UniformRand(),
                              fLo.z() + (fHi.z() - fLo.z())*G4UniformRand());
        if (fSolid->Inside(p) == kInside) ++accepted;
    }
    G4cout << "PostAct generator: volume " << name << " (" << fSolid->GetEntityType()
           << "), bounding box " << fLo/mm << " - " << fHi/mm << " mm, acceptance "
           << 100.*accepted/kTestPoints << " %" << G4endl;
}
