#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

// Mandatory user class that defines the detector used in the
// simulation, its geometry and its materials.
// Derived from the G4VUserDetectorConstruction initialisation 
// abstract base class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction();		// Constructor
		~DetectorConstruction();	// Destructor

	private:
        // Method to construct the detector
        G4VPhysicalVolume* Construct();
        // Method to define the materials
        void DefineMaterials();

	public:
		// Get methods
		// Method to get the world physical volume
	    const G4VPhysicalVolume* GetWorld()     {return experimentalHall_phys;};
			
		// Geometry update
		void UpdateGeometry();


	private:
		// Physical volumes
		G4VPhysicalVolume* experimentalHall_phys;

};

#endif

