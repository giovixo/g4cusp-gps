#include "DetectorConstruction.hh"

#include "SensitiveDetector.hh"
//#include "SDDSensitiveDetector.hh"

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

//#include "ConfigFile.hh"

// Constructor 
DetectorConstruction::DetectorConstruction() 
:	experimentalHall_log(0),
    sdd00_log(0),
    sdd01_log(0),
    sdd02_log(0),
    sdd03_log(0),
    sdd04_log(0),
    sdd05_log(0),
    sdd06_log(0),
    sdd07_log(0),
    sdd08_log(0),
    sdd09_log(0),
    sdd10_log(0),
    sdd11_log(0),
    sdd12_log(0),
    sdd13_log(0),
    sdd14_log(0),
    sdd15_log(0),
    sdd16_log(0),
    sdd17_log(0),
    sdd18_log(0),
    sdd19_log(0),
    sdd20_log(0),
    sdd21_log(0),
    sdd22_log(0),
    sdd23_log(0),
    sdd24_log(0),
    sdd25_log(0),
    sdd26_log(0),
    sdd27_log(0),
    sdd28_log(0),
    sdd29_log(0),
    sdd30_log(0),
    sdd31_log(0),

    scint00_log(0),
    scint01_log(0),
    scint02_log(0),
    scint03_log(0),
    scint04_log(0),
    scint05_log(0),
    scint06_log(0),
    scint07_log(0),
    scint08_log(0),
    scint09_log(0),
    scint10_log(0),
    scint11_log(0),
    scint12_log(0),
    scint13_log(0),
    scint14_log(0),
    scint15_log(0),
    scint16_log(0),
    scint17_log(0),
    scint18_log(0),
    scint19_log(0),
    scint20_log(0),
    scint21_log(0),
    scint22_log(0),
    scint23_log(0),
    scint24_log(0),
    scint25_log(0),
    scint26_log(0),
    scint27_log(0),
    scint28_log(0),
    scint29_log(0),
    scint30_log(0),
    scint31_log(0),
    scint32_log(0),
    scint33_log(0),
    scint34_log(0),
    scint35_log(0),
    scint36_log(0),
    scint37_log(0),
    scint38_log(0),
    scint39_log(0),
    scint40_log(0),
    scint41_log(0),
    scint42_log(0),
    scint43_log(0),
    scint44_log(0),
    scint45_log(0),
    scint46_log(0),
    scint47_log(0),
    scint48_log(0),
    scint49_log(0),
    scint50_log(0),
    scint51_log(0),
    scint52_log(0),
    scint53_log(0),
    scint54_log(0),
    scint55_log(0),
    scint56_log(0),
    scint57_log(0),
    scint58_log(0),
    scint59_log(0),
    scint60_log(0),
    scint61_log(0),
    scint62_log(0),
    scint63_log(0),
    experimentalHall_phys(0)
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
    G4Element*  Fe = new G4Element("Iron"      , "Fe", z = 26., a =  55.84*g/mole);
    G4Element*  Co = new G4Element("Cobalt"    , "Co", z = 27., a =  58.93*g/mole);
    G4Element*  Ni = new G4Element("Nickel"    , "Ni", z = 28., a =  58.69*g/mole);
//	G4Element*  Cu = new G4Element("Copper"    , "Cu", z = 29., a =  63.55*g/mole);
    G4Element*  Ga = new G4Element("Gallium"   , "Ga", z = 31., a =  69.72*g/mole);
//  G4Element*  Ge = new G4Element("Germanium" , "Ge", z = 32., a =  72.63*g/mole);
//  G4Element*  As = new G4Element("Arsenic"   , "As", z = 33., a =  74.92*g/mole);
    G4Element*  Br = new G4Element("Bromine"   , "Br", z = 35., a =  79.90*g/mole);
//  G4Element*  Mo = new G4Element("Molybdenum", "Mo", z = 42., a =  95.94*g/mole);
//  G4Element*  I  = new G4Element("Iodine"    , "I" , z = 53., a = 126.90*g/mole);
//  G4Element*  Cs = new G4Element("Caesium"   , "Cs", z = 55., a = 132.90*g/mole);
    G4Element*  La = new G4Element("Lanthanium", "La", z = 57., a = 138.90*g/mole);
    G4Element*  Gd = new G4Element("Gadolinium", "Gd", z = 64., a = 157.25*g/mole);
//  G4Element*  Ta = new G4Element("Tantalum"  , "Ta", z = 73., a = 180.94*g/mole);
    G4Element*  W  = new G4Element("Tungsten"  , "W" , z = 74., a = 183.84*g/mole);
//  G4Element*  Pb = new G4Element("Lead"      , "Pb", z = 82., a = 207.20*g/mole);
//  G4Element*  Bi = new G4Element("Bismuth"   , "Bi", z = 83., a = 208.98*g/mole);

	// Materials 
	// Vacuum
	G4Material* Vacuum = new G4Material("Vacuum", density = 1.e-25*g/cm3, nel = 1);
	Vacuum -> AddElement(H, 100*perCent);

    // NIST
    G4NistManager* man = G4NistManager::Instance();
    G4Material* G4_CESIUM_IODIDE = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
    G4Material* G4_BGO = man->FindOrBuildMaterial("G4_BGO");
    G4Material* G4_Si = man->FindOrBuildMaterial("G4_Si");
    G4Material* G4_Ti = man->FindOrBuildMaterial("G4_Ti");
    G4Material* G4_Al = man->FindOrBuildMaterial("G4_Al");
    G4Material* G4_AIR = man->FindOrBuildMaterial("G4_AIR");
    G4Material* G4_STAINLESS_STEEL = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* G4_BRASS = man->FindOrBuildMaterial("G4_BRASS");
    G4Material* G4_Cu = man->FindOrBuildMaterial("G4_Cu");
    G4Material* G4_Mo = man->FindOrBuildMaterial("G4_Mo");
    G4Material* G4_W = man->FindOrBuildMaterial("G4_W");
    G4Material* G4_Ta = man->FindOrBuildMaterial("G4_Ta");
    G4Material* G4_Pb = man->FindOrBuildMaterial("G4_Pb");
    G4Material* G4_SILICON_DIOXIDE = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* G4_KAPTON = man->FindOrBuildMaterial("G4_KAPTON");
    G4Material* G4_Galactic = man->FindOrBuildMaterial("G4_Galactic");
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
    
    // Kovar
    G4Material* Kovar = new G4Material("Kovar", density = 8.0 *g/cm3, nel = 3);
    Kovar -> AddElement(Fe, 54*perCent);
    Kovar -> AddElement(Ni, 29*perCent);
    Kovar -> AddElement(Co, 17*perCent);


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

    // W with 90 % density (=0.9*19.3)
    G4Material* W90 = new G4Material("W90", density = 17.37*g/cm3, nel = 1);
    W90 -> AddElement(W, 100*perCent);

    // Al with 10 % density (=0.1*2.699)
    G4Material* Al10 = new G4Material("Al10", density = 0.2699*g/cm3, nel = 1);
    Al10 -> AddElement(Al, 100*perCent);
    
    // Silicone (Dowsil 93-500, Polydimethylsiloxane C2H6OSi)
    G4Material* Silicone = new G4Material("Silicone", density = 1.08 *g/cm3, nel = 4);
    Silicone -> AddElement(C,  natoms=2);
    Silicone -> AddElement(H,  natoms=6);
    Silicone -> AddElement(O,  natoms=1);
    Silicone -> AddElement(Si, natoms=1);
    
    // Silicon (for SDD)
    // WARNING: in the GDML file the SDD overall thickness is 400 um instead of 450 um!
    // So we scale the 2.33*g/cm3 Si  density by a 450/400 factor, to ensure a correct efficiency
    G4Material* Si_SDD = new G4Material("Si_SDD", density = 2.62 *g/cm3, ncomponents=1);
    Si_SDD -> AddMaterial(G4_Si, fractionmass=100*perCent);


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
    
    sdd00_log = parser.GetVolume("Absorber_001");
    sdd01_log = parser.GetVolume("Absorber_002");
    sdd02_log = parser.GetVolume("Absorber_003");
    sdd03_log = parser.GetVolume("Absorber_004");
    sdd04_log = parser.GetVolume("Absorber_005");
    sdd05_log = parser.GetVolume("Absorber_006");
    sdd06_log = parser.GetVolume("Absorber_007");
    sdd07_log = parser.GetVolume("Absorber_008");
    sdd08_log = parser.GetVolume("Absorber_009");
    sdd09_log = parser.GetVolume("Absorber_010");
    sdd10_log = parser.GetVolume("Absorber_011");
    sdd11_log = parser.GetVolume("Absorber_012");
    sdd12_log = parser.GetVolume("Absorber_013");
    sdd13_log = parser.GetVolume("Absorber_014");
    sdd14_log = parser.GetVolume("Absorber_015");
    sdd15_log = parser.GetVolume("Absorber_016");
    sdd16_log = parser.GetVolume("Absorber_017");
    sdd17_log = parser.GetVolume("Absorber_018");
    sdd18_log = parser.GetVolume("Absorber_019");
    sdd19_log = parser.GetVolume("Absorber_020");
    sdd20_log = parser.GetVolume("Absorber_021");
    sdd21_log = parser.GetVolume("Absorber_022");
    sdd22_log = parser.GetVolume("Absorber_023");
    sdd23_log = parser.GetVolume("Absorber_024");
    sdd24_log = parser.GetVolume("Absorber_025");
    sdd25_log = parser.GetVolume("Absorber_026");
    sdd26_log = parser.GetVolume("Absorber_027");
    sdd27_log = parser.GetVolume("Absorber_028");
    sdd28_log = parser.GetVolume("Absorber_029");
    sdd29_log = parser.GetVolume("Absorber_030");
    sdd30_log = parser.GetVolume("Absorber_031");
    sdd31_log = parser.GetVolume("Absorber_032");

    scint00_log = parser.GetVolume("Scatterer_001");
    scint01_log = parser.GetVolume("Scatterer_002");
    scint02_log = parser.GetVolume("Scatterer_003");
    scint03_log = parser.GetVolume("Scatterer_004");
    scint04_log = parser.GetVolume("Scatterer_005");
    scint05_log = parser.GetVolume("Scatterer_006");
    scint06_log = parser.GetVolume("Scatterer_007");
    scint07_log = parser.GetVolume("Scatterer_008");
    scint08_log = parser.GetVolume("Scatterer_009");
    scint09_log = parser.GetVolume("Scatterer_010");
    scint10_log = parser.GetVolume("Scatterer_011");
    scint11_log = parser.GetVolume("Scatterer_012");
    scint12_log = parser.GetVolume("Scatterer_013");
    scint13_log = parser.GetVolume("Scatterer_014");
    scint14_log = parser.GetVolume("Scatterer_015");
    scint15_log = parser.GetVolume("Scatterer_016");
    scint16_log = parser.GetVolume("Scatterer_017");
    scint17_log = parser.GetVolume("Scatterer_018");
    scint18_log = parser.GetVolume("Scatterer_019");
    scint19_log = parser.GetVolume("Scatterer_020");
    scint20_log = parser.GetVolume("Scatterer_021");
    scint21_log = parser.GetVolume("Scatterer_022");
    scint22_log = parser.GetVolume("Scatterer_023");
    scint23_log = parser.GetVolume("Scatterer_024");
    scint24_log = parser.GetVolume("Scatterer_025");
    scint25_log = parser.GetVolume("Scatterer_026");
    scint26_log = parser.GetVolume("Scatterer_027");
    scint27_log = parser.GetVolume("Scatterer_028");
    scint28_log = parser.GetVolume("Scatterer_029");
    scint29_log = parser.GetVolume("Scatterer_030");
    scint30_log = parser.GetVolume("Scatterer_031");
    scint31_log = parser.GetVolume("Scatterer_032");
    scint32_log = parser.GetVolume("Scatterer_033");
    scint33_log = parser.GetVolume("Scatterer_034");
    scint34_log = parser.GetVolume("Scatterer_035");
    scint35_log = parser.GetVolume("Scatterer_036");
    scint36_log = parser.GetVolume("Scatterer_037");
    scint37_log = parser.GetVolume("Scatterer_038");
    scint38_log = parser.GetVolume("Scatterer_039");
    scint39_log = parser.GetVolume("Scatterer_040");
    scint40_log = parser.GetVolume("Scatterer_041");
    scint41_log = parser.GetVolume("Scatterer_042");
    scint42_log = parser.GetVolume("Scatterer_043");
    scint43_log = parser.GetVolume("Scatterer_044");
    scint44_log = parser.GetVolume("Scatterer_045");
    scint45_log = parser.GetVolume("Scatterer_046");
    scint46_log = parser.GetVolume("Scatterer_047");
    scint47_log = parser.GetVolume("Scatterer_048");
    scint48_log = parser.GetVolume("Scatterer_049");
    scint49_log = parser.GetVolume("Scatterer_050");
    scint50_log = parser.GetVolume("Scatterer_051");
    scint51_log = parser.GetVolume("Scatterer_052");
    scint52_log = parser.GetVolume("Scatterer_053");
    scint53_log = parser.GetVolume("Scatterer_054");
    scint54_log = parser.GetVolume("Scatterer_055");
    scint55_log = parser.GetVolume("Scatterer_056");
    scint56_log = parser.GetVolume("Scatterer_057");
    scint57_log = parser.GetVolume("Scatterer_058");
    scint58_log = parser.GetVolume("Scatterer_059");
    scint59_log = parser.GetVolume("Scatterer_060");
    scint60_log = parser.GetVolume("Scatterer_061");
    scint61_log = parser.GetVolume("Scatterer_062");
    scint62_log = parser.GetVolume("Scatterer_063");
    scint63_log = parser.GetVolume("Scatterer_064");

    
	// The function must return the physical volume of the world
	return experimentalHall_phys;
}


void DetectorConstruction::ConstructSDandField()
{
    // Sensitive volume
    // So far, the detector created above is not yet a real detector: it is just
    // a geometrical object (with some attributes) placed within the world volume.
    // To make the volume a detector, which can record e.g. hits, one must
    // define a sensitive volume associated with it.
    // For this purpose, a SensitiveDetector object is instantiated.
    // SD for the Scint
    

    auto sdman = G4SDManager::GetSDMpointer(); // Mandatory since Geant v. 4.10.03

    // Instantiation of the scintillator sensitive detector and readout geometry
    SensitiveDetector* scint_SD  = new SensitiveDetector("SCI");


    sdman->AddNewDetector(scint_SD); // Mandatory since Geant v. 4.10.03

    SetSensitiveDetector(sdd00_log, scint_SD);
    SetSensitiveDetector(sdd01_log, scint_SD);
    SetSensitiveDetector(sdd02_log, scint_SD);
    SetSensitiveDetector(sdd03_log, scint_SD);
    SetSensitiveDetector(sdd04_log, scint_SD);
    SetSensitiveDetector(sdd05_log, scint_SD);
    SetSensitiveDetector(sdd06_log, scint_SD);
    SetSensitiveDetector(sdd07_log, scint_SD);
    SetSensitiveDetector(sdd08_log, scint_SD);
    SetSensitiveDetector(sdd09_log, scint_SD);
    SetSensitiveDetector(sdd10_log, scint_SD);
    SetSensitiveDetector(sdd11_log, scint_SD);
    SetSensitiveDetector(sdd12_log, scint_SD);
    SetSensitiveDetector(sdd13_log, scint_SD);
    SetSensitiveDetector(sdd14_log, scint_SD);
    SetSensitiveDetector(sdd15_log, scint_SD);
    SetSensitiveDetector(sdd16_log, scint_SD);
    SetSensitiveDetector(sdd17_log, scint_SD);
    SetSensitiveDetector(sdd18_log, scint_SD);
    SetSensitiveDetector(sdd19_log, scint_SD);
    SetSensitiveDetector(sdd20_log, scint_SD);
    SetSensitiveDetector(sdd21_log, scint_SD);
    SetSensitiveDetector(sdd22_log, scint_SD);
    SetSensitiveDetector(sdd23_log, scint_SD);
    SetSensitiveDetector(sdd24_log, scint_SD);
    SetSensitiveDetector(sdd25_log, scint_SD);
    SetSensitiveDetector(sdd26_log, scint_SD);
    SetSensitiveDetector(sdd27_log, scint_SD);
    SetSensitiveDetector(sdd28_log, scint_SD);
    SetSensitiveDetector(sdd29_log, scint_SD);
    SetSensitiveDetector(sdd30_log, scint_SD);
    SetSensitiveDetector(sdd31_log, scint_SD);

    SetSensitiveDetector(scint00_log, scint_SD);
    SetSensitiveDetector(scint01_log, scint_SD);
    SetSensitiveDetector(scint02_log, scint_SD);
    SetSensitiveDetector(scint03_log, scint_SD);
    SetSensitiveDetector(scint04_log, scint_SD);
    SetSensitiveDetector(scint05_log, scint_SD);
    SetSensitiveDetector(scint06_log, scint_SD);
    SetSensitiveDetector(scint07_log, scint_SD);
    SetSensitiveDetector(scint08_log, scint_SD);
    SetSensitiveDetector(scint09_log, scint_SD);
    SetSensitiveDetector(scint10_log, scint_SD);
    SetSensitiveDetector(scint11_log, scint_SD);
    SetSensitiveDetector(scint12_log, scint_SD);
    SetSensitiveDetector(scint13_log, scint_SD);
    SetSensitiveDetector(scint14_log, scint_SD);
    SetSensitiveDetector(scint15_log, scint_SD);
    SetSensitiveDetector(scint16_log, scint_SD);
    SetSensitiveDetector(scint17_log, scint_SD);
    SetSensitiveDetector(scint18_log, scint_SD);
    SetSensitiveDetector(scint19_log, scint_SD);
    SetSensitiveDetector(scint20_log, scint_SD);
    SetSensitiveDetector(scint21_log, scint_SD);
    SetSensitiveDetector(scint22_log, scint_SD);
    SetSensitiveDetector(scint23_log, scint_SD);
    SetSensitiveDetector(scint24_log, scint_SD);
    SetSensitiveDetector(scint25_log, scint_SD);
    SetSensitiveDetector(scint26_log, scint_SD);
    SetSensitiveDetector(scint27_log, scint_SD);
    SetSensitiveDetector(scint28_log, scint_SD);
    SetSensitiveDetector(scint29_log, scint_SD);
    SetSensitiveDetector(scint30_log, scint_SD);
    SetSensitiveDetector(scint31_log, scint_SD);
    SetSensitiveDetector(scint32_log, scint_SD);
    SetSensitiveDetector(scint33_log, scint_SD);
    SetSensitiveDetector(scint34_log, scint_SD);
    SetSensitiveDetector(scint35_log, scint_SD);
    SetSensitiveDetector(scint36_log, scint_SD);
    SetSensitiveDetector(scint37_log, scint_SD);
    SetSensitiveDetector(scint38_log, scint_SD);
    SetSensitiveDetector(scint39_log, scint_SD);
    SetSensitiveDetector(scint40_log, scint_SD);
    SetSensitiveDetector(scint41_log, scint_SD);
    SetSensitiveDetector(scint42_log, scint_SD);
    SetSensitiveDetector(scint43_log, scint_SD);
    SetSensitiveDetector(scint44_log, scint_SD);
    SetSensitiveDetector(scint45_log, scint_SD);
    SetSensitiveDetector(scint46_log, scint_SD);
    SetSensitiveDetector(scint47_log, scint_SD);
    SetSensitiveDetector(scint48_log, scint_SD);
    SetSensitiveDetector(scint49_log, scint_SD);
    SetSensitiveDetector(scint50_log, scint_SD);
    SetSensitiveDetector(scint51_log, scint_SD);
    SetSensitiveDetector(scint52_log, scint_SD);
    SetSensitiveDetector(scint53_log, scint_SD);
    SetSensitiveDetector(scint54_log, scint_SD);
    SetSensitiveDetector(scint55_log, scint_SD);
    SetSensitiveDetector(scint56_log, scint_SD);
    SetSensitiveDetector(scint57_log, scint_SD);
    SetSensitiveDetector(scint58_log, scint_SD);
    SetSensitiveDetector(scint59_log, scint_SD);
    SetSensitiveDetector(scint60_log, scint_SD);
    SetSensitiveDetector(scint61_log, scint_SD);
    SetSensitiveDetector(scint62_log, scint_SD);
    SetSensitiveDetector(scint63_log, scint_SD);    
}




void DetectorConstruction::UpdateGeometry()
{  
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


