#ifndef PostActActionInitialization_h
#define PostActActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

// User actions of cusp-postactivation: decay of the activated nuclides in the
// mass model, one (volume, isotope) pair per run.
// Work in progress: only a placeholder primary generator is registered.

class PostActActionInitialization : public G4VUserActionInitialization
{
public:
    PostActActionInitialization();
    virtual ~PostActActionInitialization();

    virtual void Build() const;
    virtual void BuildForMaster() const;
};

#endif
