#include "PostActStackingAction.hh"
#include "PostActConfig.hh"

#include "G4Ions.hh"
#include "G4Track.hh"

#include <cmath>


PostActStackingAction::PostActStackingAction()
: G4UserStackingAction()
{;}


PostActStackingAction::~PostActStackingAction()
{;}


G4bool PostActStackingAction::IsLongLived(const G4ParticleDefinition* particle)
{
    if (dynamic_cast<const G4Ions*>(particle) == nullptr) return false;
    if (particle->GetPDGStable()) return false;
    const G4double lifetime = particle->GetPDGLifeTime();
    return lifetime > 0. && lifetime*std::log(2.) >= PostActConfig::MinHalfLife();
}


G4ClassificationOfNewTrack PostActStackingAction::ClassifyNewTrack(const G4Track* track)
{
    // Daughter nuclide created at rest and not to decay here
    if (track->GetParentID() > 0 && track->GetKineticEnergy() == 0. &&
        IsLongLived(track->GetDefinition()))
        return fKill;
    return fUrgent;
}
