#ifndef DETECTORHIT_HH
#define DETECTORHIT_HH

#include <vector>
#include <utility>
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "tls.hh"

// General purpose:
//  Represents a hit and contains relevant information about the hit.
//  (The user implementation should provide functionalities which allow to 
//  assign and retrieve hit information).
//  Hit objects can be created in the sensitive detector, and may be used
//  at the end of events to accumulate hit information.
//
// Purpose of the class in this implementation:
//  Represents a hit in the detector, where the relevant hit information is 
//  the energy deposited by the particle and its location


// G4VHit is the abstract base class for creating hit objects.

class DetectorHit : public G4VHit {
	public:
		DetectorHit();							// Default constructor
		DetectorHit(G4int xp);                  // Alternate constructor
		virtual ~DetectorHit();					// Destructor

		
		// Assignment and comparison operators:  
		const DetectorHit& operator= (const DetectorHit& right);
		int operator==(const DetectorHit& right) const;

		
		// The hit class has user-defined new and delete operators to speed up
		// the generation and deletion of hits objects:
		inline void* operator new(size_t);
		inline void  operator delete(void* hit);

		
		// The G4VHit provides two methods, which can be overloaded by the user to
		// visualize hits or to print information about hits. 
		// Here, these methods are not used (dummy implementation):
		virtual void Draw() { }
		virtual void Print() { }

		
		// The current hit object can be used to keep track of:
		// a. The  energy deposit
		// b. The ID of the sector of deposit
        // c. The coordinates of the voxel in which the energy deposit happened
		
		// The following methods allow to set, accumulate and get the energy deposit:
		inline void SetEnergyDeposit(G4double energy) {energyDeposit = energy;}
		inline void AddEnergyDeposit(G4double energy) {energyDeposit += energy;}
		inline G4double GetEnergyDeposit() 			  {return energyDeposit;}
		
		// The following methods allow to set and get the ID of the scintillator bar:
		inline void SetScintID(G4int id) 			  {scintID = id;}
		inline G4int GetScintID() 					  {return scintID;}

		// The following methods allow to set and get the position (x,y,z) ID of the scintillator bar:
		inline void SetScint_x(G4double x) 			  {scint_x = x;}
		inline void SetScint_y(G4double y) 			  {scint_y = y;}
		inline void SetScint_z(G4double z) 			  {scint_z = z;}
		inline G4double GetScint_x()   {return scint_x;} 
		inline G4double GetScint_y()   {return scint_y;}
		inline G4double GetScint_z()   {return scint_z;}
	
		// The following methods allow to set and get the time:
		inline void SetTime(G4double t) 		      {globaltime = t;}
		inline G4int GetTime() 					      {return globaltime;}

	private:
		G4double energyDeposit;
		G4int scintID;
		G4double globaltime;
		G4double scint_x;
		G4double scint_y;
		G4double scint_z;
	
};


typedef G4THitsCollection<DetectorHit> DetectorHitsCollection;
extern G4ThreadLocal G4Allocator<DetectorHit>* DetectorHitAllocator;


inline void* DetectorHit::operator new(size_t) {
    // Implementation of the new operator
    if(!DetectorHitAllocator) DetectorHitAllocator = new G4Allocator<DetectorHit>;
    return (void*) DetectorHitAllocator->MallocSingle();
}


inline void DetectorHit::operator delete(void* hit) {
    // Implementation of the delete operator
    DetectorHitAllocator->FreeSingle((DetectorHit*) hit);
}

#endif // DETECTORHIT_HH
