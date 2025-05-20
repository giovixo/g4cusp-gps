#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include "DetectorHit.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


// G4VSensitiveDetector is the abstract base class for sensitive detectors.

class SensitiveDetector : public G4VSensitiveDetector 
{
	public:
        SensitiveDetector(G4String SDuniqueName);
        // Constructor w/name of SD and array dim.
        virtual ~SensitiveDetector();				    			// Destructor

		// The G4VSensitiveDetector class has two methods, which may be overloaded
		// by the user: Initialize and EndOfEvent. These two methods are invoked
		// at the beginning and the end of each event.
		virtual void Initialize(G4HCofThisEvent*);
		virtual void EndOfEvent(G4HCofThisEvent*);

	protected:
		// The principal purpose of a sensitive detector is to construct hit objects
		// using information from particle steps. 
		// These hits should be created in the ProcessHits function, which is called
		// each time a particle performs a step inside a logical volume which has a 
		// sensitive detector attached.
		// Two arguments are passed to the ProcessHits method: 
		// 1) The first one is the current G4Step object, which can be used to 
		//    retrieve information about the step like energy loss, physics 
		//    processes,... 
		// 2) The second argument is a G4TouchableHistory object for the read-out 
		//    geometry (if a read-out geometry was assigned to the sensitive 
		//    detector), which can be used e.g. to determine read-out channels,... 
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* hist);

	private:
		DetectorHitsCollection* hitsCollection;
		G4int hitsCollectionID;
    
        const int numberOfCellsX;
        const int numberOfCellsY;
        int **CellID;
        void initArray();
};

#endif // SENSITIVEDETECTOR_HH
