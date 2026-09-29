#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;

#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"

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
        // Method to read the values of parameters from file
        void DefineParameters();
        // My messanger
        G4GenericMessenger *fMessenger;

	public:
		// Get methods
		// Method to get the world physical volume
	    const G4VPhysicalVolume* GetWorld()     {return experimentalHall_phys;};
			
		// Geometry update
		void UpdateGeometry();


	private:
        // Outout (root) file name
        G4String file_name;


		// Logical volumes
    

        

		// Physical volumes
		G4VPhysicalVolume* experimentalHall_phys;

        G4Material* mli1Material;
        G4Material* mli2Material;
        G4Material* mli3Material;

        G4Material* solarPanel1Material;
        G4Material* solarPanel2Material;
        G4Material* busMaterial;

};

#endif

