#ifndef PostActActionInitialization_h
#define PostActActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

// User actions of cusp-postactivation: decay of the activated nuclides in the
// mass model, one (volume, isotope) pair per run.
// Registers the generator, run, event, stepping and stacking actions.

class PostActActionInitialization : public G4VUserActionInitialization
{
public:
    PostActActionInitialization();
    virtual ~PostActActionInitialization();

    virtual void Build() const;
    virtual void BuildForMaster() const;
};

#endif
