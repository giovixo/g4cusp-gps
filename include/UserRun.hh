#ifndef USERRUN_HH
#define USERRUN_HH

#include "G4Run.hh"
#include "G4Event.hh"
#include "globals.hh"

#include "G4GenericMessenger.hh"

// To collect data over worker threads at the end of an event loop, derived class of G4Run should be used. 
// The base class defines the following virtual methods:
// void RecordEvent(const G4Event*)
//      Method to be overwritten by the user for recording events in this (thread-local) run.
//      Invoke G4Run base-class method for recording data member in the base class.
//      If the user in the past has implemented an access to hits or scores in his/her UserEventAction::EndOfEventAction(), 
//      (s)he can simply move that code segment to this RecordEvent() method. See the following sample code.
// void Merge(const G4Run*)
//      Method to be overwritten by the user for merging local Run object to the global Run object.
//      Invoke G4Run base-class method to merge data member in the base class.


class UserRun : public G4Run
{
    public:
        UserRun();
        virtual ~UserRun();
        void RecordEvent(const G4Event*);
        void Merge(const G4Run*);
        
    private:
        G4int SDD_hitsCollectionIndex;
        G4int SCI_hitsCollectionIndex;

        // Output (root) file name
        G4String file_name;
        // My messanger
        //G4GenericMessenger *fMessenger;
};

#endif // USERRUN_HH
