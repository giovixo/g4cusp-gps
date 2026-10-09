#ifndef PostActPrimaryGeneratorAction_h
#define PostActPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "globals.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class G4VPhysicalVolume;
class G4VSolid;
class G4Event;

// Primary generator of cusp-postactivation.
// Each event is one nucleus of the isotope set by /postact/isotope, at rest, at a point
// uniformly distributed inside the physical volume set by /postact/volume.
// The isotope and the volume are re-read at every event and the derived data are cached
// per thread, so they can be changed between runs.

class PostActPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PostActPrimaryGeneratorAction();
    virtual ~PostActPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* anEvent);

private:
    void UpdateIon();                            // (re)creates the ion if the isotope changed
    void UpdateVolume();                         // (re)reads the volume geometry if the volume changed
    G4bool Sample(G4ThreeVector& local) const;   // one uniform point in the solid (local frame)

    // Placement chain from mother down to target; transform = local(target) -> frame of mother
    static G4bool FindPath(const G4VPhysicalVolume* mother, const G4VPhysicalVolume* target,
                           G4Transform3D& transform);

    G4ParticleGun* fParticleGun;

    // Ion cache
    G4String fCurrentIsotope;
    G4ParticleDefinition* fIon;

    // Volume cache
    G4String fCurrentVolume;
    const G4VSolid* fSolid;
    G4Transform3D fTransform;                    // local frame of the volume -> global frame
    G4ThreeVector fLo, fHi;                      // bounding box of the solid (local frame)
};

#endif
