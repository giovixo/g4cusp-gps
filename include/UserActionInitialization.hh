#ifndef UserActionInitialization_h
#define UserActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class UserActionInitialization : public G4VUserActionInitialization
{
public:
    UserActionInitialization();	            // Constructor
    virtual ~UserActionInitialization();    // Destructor
    
    virtual void Build() const;
    virtual void BuildForMaster() const;
};

#endif
