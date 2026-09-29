#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4RegionStore.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Trap.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"
#include "globals.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4GDMLParser.hh"

#include <fstream>
#include <vector>
#include <string>


// Constructor 
DetectorConstruction::DetectorConstruction() 
:   experimentalHall_phys(0)
{
    // User messanger
    // fMessenger = new G4GenericMessenger(this, "/parameters/", "Output file name");
    // fMessenger->DeclareProperty("filename", file_name, "Name of the output file");

	// Define the materials
	DefineMaterials();

	// Define the parameters
	DefineParameters();
		

	
;}



// Destructor
DetectorConstruction::~DetectorConstruction() { }



// Definition of the parameters
void DetectorConstruction::DefineParameters()
{
}


// Definition of the materials
void DetectorConstruction::DefineMaterials()
{
	G4double a;			// Atomic mass
	G4double z;			// Atomic number
	G4double density;	// Density
	G4int nel;			// Number of elements in a compound
	G4int natoms;       // Number of atoms in a compound
	G4double fractionmass;
    G4int ncomponents;

	// Elements
	G4Element*  H  = new G4Element("Hydrogen"  , "H" , z = 1. , a =  1.008*g/mole);
//	G4Element*  He = new G4Element("Helium"    , "He", z = 2. , a =  4.003*g/mole);
    G4Element*  C  = new G4Element("Carbon"    , "C" , z = 6. , a =  12.01*g/mole);
    G4Element*  N  = new G4Element("Nitrogen"  , "N" , z = 7. , a =  14.01*g/mole);
    G4Element*  O  = new G4Element("Oxygen"    , "O" , z = 8. , a =  16.00*g/mole);
//  G4Element*  Na = new G4Element("Sodium"    , "Na", z = 11., a =  22.99*g/mole);
	G4Element*  Al = new G4Element("Aluminium" , "Al", z = 13., a =  26.98*g/mole);
	G4Element*  Si = new G4Element("Silicon"   , "Si", z = 14., a =  28.08*g/mole);
//  G4Element*  K  = new G4Element("Potassium" , "K" , z = 19., a =  39.10*g/mole);
//  G4Element*  Ti = new G4Element("Titanium"  , "Ti", z = 22., a =  47.87*g/mole);
//	G4Element*  Cu = new G4Element("Copper"    , "Cu", z = 29., a =  63.55*g/mole);
    G4Element*  Ga = new G4Element("Gallium"   , "Ga", z = 31., a =  69.72*g/mole);
//  G4Element*  Ge = new G4Element("Germanium" , "Ge", z = 32., a =  72.63*g/mole);
//  G4Element*  As = new G4Element("Arsenic"   , "As", z = 33., a =  74.92*g/mole);
//  G4Element*  Mo = new G4Element("Molybdenum", "Mo", z = 42., a =  95.94*g/mole);
//  G4Element*  I  = new G4Element("Iodine"    , "I" , z = 53., a = 126.90*g/mole);
//  G4Element*  Cs = new G4Element("Caesium"   , "Cs", z = 55., a = 132.90*g/mole);
    G4Element*  Gd = new G4Element("Gadolinium", "Gd", z = 64., a = 157.25*g/mole);
//  G4Element*  Ta = new G4Element("Tantalum"  , "Ta", z = 73., a = 180.94*g/mole);
//  G4Element*  Pb = new G4Element("Lead"      , "Pb", z = 82., a = 207.20*g/mole);
//  G4Element*  Bi = new G4Element("Bismuth"   , "Bi", z = 83., a = 208.98*g/mole);

	// Materials 
	// Vacuum
	G4Material* Vacuum = new G4Material("Vacuum", density = 1.e-25*g/cm3, nel = 1);
	Vacuum -> AddElement(H, 100*perCent);

    // NIST
    G4NistManager* man = G4NistManager::Instance();
    G4Material* G4_Al = man->FindOrBuildMaterial("G4_Al");
    G4Material* G4_SILICON_DIOXIDE = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* G4_KAPTON = man->FindOrBuildMaterial("G4_KAPTON");
    G4Material* G4_POLYPROPYLENE = man->FindOrBuildMaterial("G4_POLYPROPYLENE");

    // Effective Aluminium Solid for Bus (6619.5 g before system margin, minus 1579.2 g payload = 5040 g)
    // Volume bus: 10 x 10 x 20 = 2000 cm3; effective density = 2.52 g/cm3
    G4Material* EffectiveAluminiumSolid_Bus = new G4Material("EffectiveAluminiumSolid_Bus",
                                                         density = 2.52*g/cm3, nel = 1);
    EffectiveAluminiumSolid_Bus -> AddElement(Al, 100*perCent);
    
    // Effective Aluminium Solid for CSAC (35 g, 17 cm3)
    G4Material* EffectiveAluminiumSolid_CSAC = new G4Material("EffectiveAluminiumSolid_CSAC",
                                                         density = 2.06*g/cm3, nel = 1);
    EffectiveAluminiumSolid_CSAC -> AddMaterial(G4_Al, 100*perCent);
    
    
    // Effective Aluminium Solid for PICO SAR 250 DC-DC converter (12 g, 1.10" x 0.80" x 0.45")
    G4Material* EffectiveAluminiumSolid_PICO = new G4Material("EffectiveAluminiumSolid_PICO",
                                                         density = 1.85*g/cm3, nel = 1);
    EffectiveAluminiumSolid_PICO -> AddMaterial(G4_Al, 100*perCent);


    // GAGG
    G4Material* GAGG = new G4Material("GAGG", density = 6.63 *g/cm3, nel = 4);
    GAGG -> AddElement(Gd, natoms=3);
    GAGG -> AddElement(Al, natoms=2);
    GAGG -> AddElement(Ga, natoms=3);
    GAGG -> AddElement(O,  natoms=12);
    

    G4Material* OpticalFilter = new G4Material("OpticalFilter", density = 0.2223*g/cm3, ncomponents=2);
    OpticalFilter -> AddMaterial(G4_KAPTON, fractionmass=63.2*perCent);
    OpticalFilter -> AddMaterial(G4_Al,  fractionmass=36.8*perCent);

    // Diglycidyl Ether of Bisphenol A (First compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_1 = new G4Material("Epoxy_1", density = 1.16*g/cm3, nel = 3);
    Epoxy_1 -> AddElement(C, natoms=19);
    Epoxy_1 -> AddElement(H, natoms=20);
    Epoxy_1 -> AddElement(O, natoms=4);
    
    // 1,4-Butanediol Diglycidyl Ether (Second compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_2 = new G4Material("Epoxy_2", density = 1.10*g/cm3, nel = 3);
    Epoxy_2 -> AddElement(C, natoms=10);
    Epoxy_2 -> AddElement(H, natoms=18);
    Epoxy_2 -> AddElement(O, natoms=4);
    
    // 1,6-Hexanediamine 2,2,4-trimetyl (Third compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_3 = new G4Material("Epoxy_3", density = 1.16*g/cm3, nel = 3);
    Epoxy_3 -> AddElement(C, natoms=9);
    Epoxy_3 -> AddElement(H, natoms=22);
    Epoxy_3 -> AddElement(N, natoms=2);
    
    // Epoxy resin Epotek 301-1
    G4Material* Epoxy_Resin = new G4Material("Epoxy_Resin", density = 1.19*g/cm3, ncomponents = 3);
    Epoxy_Resin -> AddMaterial(Epoxy_1, fractionmass=56*perCent);
    Epoxy_Resin -> AddMaterial(Epoxy_2, fractionmass=24*perCent);
    Epoxy_Resin -> AddMaterial(Epoxy_3, fractionmass=20*perCent);
    
    // FR4 PCB material
    G4Material* FR4 = new G4Material("FR4", density = 1.8*g/cm3, ncomponents=2);
    FR4 -> AddMaterial(G4_SILICON_DIOXIDE, fractionmass=60*perCent);
    FR4 -> AddMaterial(Epoxy_Resin,  fractionmass=40*perCent);

    // Al with 10 % density (=0.1*2.699)
    G4Material* Al10 = new G4Material("Al10", density = 0.2699*g/cm3, nel = 1);
    Al10 -> AddElement(Al, 100*perCent);
    
    // Silicone (Dowsil 93-500, Polydimethylsiloxane C2H6OSi)
    G4Material* Silicone = new G4Material("Silicone", density = 1.08 *g/cm3, nel = 4);
    Silicone -> AddElement(C,  natoms=2);
    Silicone -> AddElement(H,  natoms=6);
    Silicone -> AddElement(O,  natoms=1);
    Silicone -> AddElement(Si, natoms=1);


    mli1Material = G4_KAPTON;
    mli2Material = G4_Al;
    mli3Material = G4_POLYPROPYLENE;

    solarPanel1Material = G4_Al;
    solarPanel2Material = G4_SILICON_DIOXIDE;
    busMaterial = EffectiveAluminiumSolid_Bus;
}



// Detector construction
G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// Clean old geometry, if any
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();
    
    G4GDMLParser parser;
    // Importing geometry
    //parser.Read("hermes_test.gdml");
    parser.Read("CUSP_GEANT4_Model_20240502.gdml");

    
    // Reads and stores in memory
    experimentalHall_phys = parser.GetWorldVolume(); // get world
    //experimentalHall_log = parser.GetVolume("worldVOL");
    experimentalHall_phys = parser.GetWorldVolume("CUSP_GEANT4_Model_20240502"); // get world 
    

    
	// The function must return the physical volume of the world
	return experimentalHall_phys;
}



void DetectorConstruction::UpdateGeometry()
{  
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


