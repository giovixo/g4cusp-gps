#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class SensitiveDetector;
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
        // Method to construct the sensitive detector
        void ConstructSDandField();
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
		G4LogicalVolume* experimentalHall_log;
    

        G4LogicalVolume* sdd00_log;
        G4LogicalVolume* sdd01_log;
        G4LogicalVolume* sdd02_log;
        G4LogicalVolume* sdd03_log;
        G4LogicalVolume* sdd04_log;
        G4LogicalVolume* sdd05_log;
        G4LogicalVolume* sdd06_log;
        G4LogicalVolume* sdd07_log;
        G4LogicalVolume* sdd08_log;
        G4LogicalVolume* sdd09_log;
        G4LogicalVolume* sdd10_log;
        G4LogicalVolume* sdd11_log;
        G4LogicalVolume* sdd12_log;
        G4LogicalVolume* sdd13_log;
        G4LogicalVolume* sdd14_log;
        G4LogicalVolume* sdd15_log;
        G4LogicalVolume* sdd16_log;
        G4LogicalVolume* sdd17_log;
        G4LogicalVolume* sdd18_log;
        G4LogicalVolume* sdd19_log;
        G4LogicalVolume* sdd20_log;
        G4LogicalVolume* sdd21_log;
        G4LogicalVolume* sdd22_log;
        G4LogicalVolume* sdd23_log;
        G4LogicalVolume* sdd24_log;
        G4LogicalVolume* sdd25_log;
        G4LogicalVolume* sdd26_log;
        G4LogicalVolume* sdd27_log;
        G4LogicalVolume* sdd28_log;
        G4LogicalVolume* sdd29_log;
        G4LogicalVolume* sdd30_log;
        G4LogicalVolume* sdd31_log;
        
        G4LogicalVolume* scint00_log;
        G4LogicalVolume* scint01_log;
        G4LogicalVolume* scint02_log;
        G4LogicalVolume* scint03_log;
        G4LogicalVolume* scint04_log;
        G4LogicalVolume* scint05_log;
        G4LogicalVolume* scint06_log;
        G4LogicalVolume* scint07_log;
        G4LogicalVolume* scint08_log;
        G4LogicalVolume* scint09_log;
        G4LogicalVolume* scint10_log;
        G4LogicalVolume* scint11_log;
        G4LogicalVolume* scint12_log;
        G4LogicalVolume* scint13_log;
        G4LogicalVolume* scint14_log;
        G4LogicalVolume* scint15_log;
        G4LogicalVolume* scint16_log;
        G4LogicalVolume* scint17_log;
        G4LogicalVolume* scint18_log;
        G4LogicalVolume* scint19_log;
        G4LogicalVolume* scint20_log;
        G4LogicalVolume* scint21_log;
        G4LogicalVolume* scint22_log;
        G4LogicalVolume* scint23_log;
        G4LogicalVolume* scint24_log;
        G4LogicalVolume* scint25_log;
        G4LogicalVolume* scint26_log;
        G4LogicalVolume* scint27_log;
        G4LogicalVolume* scint28_log;
        G4LogicalVolume* scint29_log;
        G4LogicalVolume* scint30_log;
        G4LogicalVolume* scint31_log;
        G4LogicalVolume* scint32_log;
        G4LogicalVolume* scint33_log;
        G4LogicalVolume* scint34_log;
        G4LogicalVolume* scint35_log;
        G4LogicalVolume* scint36_log;
        G4LogicalVolume* scint37_log;
        G4LogicalVolume* scint38_log;
        G4LogicalVolume* scint39_log;
        G4LogicalVolume* scint40_log;
        G4LogicalVolume* scint41_log;
        G4LogicalVolume* scint42_log;
        G4LogicalVolume* scint43_log;
        G4LogicalVolume* scint44_log;
        G4LogicalVolume* scint45_log;
        G4LogicalVolume* scint46_log;
        G4LogicalVolume* scint47_log;
        G4LogicalVolume* scint48_log;
        G4LogicalVolume* scint49_log;
        G4LogicalVolume* scint50_log;
        G4LogicalVolume* scint51_log;
        G4LogicalVolume* scint52_log;
        G4LogicalVolume* scint53_log;
        G4LogicalVolume* scint54_log;
        G4LogicalVolume* scint55_log;
        G4LogicalVolume* scint56_log;
        G4LogicalVolume* scint57_log;
        G4LogicalVolume* scint58_log;
        G4LogicalVolume* scint59_log;
        G4LogicalVolume* scint60_log;
        G4LogicalVolume* scint61_log;
        G4LogicalVolume* scint62_log;
        G4LogicalVolume* scint63_log;       

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

