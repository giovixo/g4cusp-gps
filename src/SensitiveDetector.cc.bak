#include "SensitiveDetector.hh"
#include "DetectorHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VSolid.hh"
#include "G4VisExtent.hh"

#include <iostream>
#include <string> 

// Constructor.
SensitiveDetector::SensitiveDetector(G4String SDuniqueName) :
G4VSensitiveDetector(SDuniqueName),
hitsCollectionID(-1),
numberOfCellsX(136),
numberOfCellsY(1)
{
    collectionName.insert(SDuniqueName);
    initArray();

}


// Destructor
SensitiveDetector::~SensitiveDetector() 
{
}

// Method for the array allocation
void SensitiveDetector::initArray()
{
    // Dynamic allocation of 2-D arrays
    // We use the sizeof() a pointer to ensure 32 and 64-bit portability
    CellID = (int **) std::calloc(numberOfCellsX, sizeof(*CellID));
    for(int iRow=0; iRow < numberOfCellsX; iRow++)
    {
        CellID[iRow] = (int *) std::calloc(numberOfCellsY, sizeof(*CellID));
    }
}



void SensitiveDetector::Initialize(G4HCofThisEvent* hitsCollEvent)
{
    // This method is called at each new event, and it gets the HCofThisEvents
    // for the event.
    hitsCollection = new DetectorHitsCollection(SensitiveDetectorName, collectionName[0]);
    
    hitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hitsCollEvent->AddHitsCollection(hitsCollectionID, hitsCollection);

    // Initialize vector of cell hits
    for(int iRow=0; iRow < numberOfCellsX; iRow++)
    {
        for(G4int iCol=0; iCol < numberOfCellsY; iCol++)
        {
            CellID[iRow][iCol] = -1;
        }
    }
}


G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory* ROhist) 
{
	// This method is called when a particle goes through a sensitive 
	// detector. It possibly creates a new hit, and add its to the collection.
	// It eventually returns a boolean value. 
	// The arguments it gets are the G4Step and the G4TouchableHistory.
	//
	// The touchable history that we get in this method contains the full information about the 
	// genealogy of the current volume in the read-out geometry. One can use the touchable to 
	// access the information on the position/rotation of a given volume.

//	if(!ROhist) return false;

    // Retrieve the energy deposited from the step
    G4double energyDeposit = step -> GetTotalEnergyDeposit();
    G4double energy_threshold = 0.05 * keV;
    if (energyDeposit < energy_threshold) return false;
    
    G4Track* track = step -> GetTrack();
    G4double globaltime = track -> GetGlobalTime();

    G4double x_center;
    G4double y_center;
    G4double z_center;
    G4VSolid *targetSolid = 0;
    G4VisExtent targetVisExtent;
    
    
    // debug code
//    G4StepPoint* point1 = step->GetPreStepPoint();
//    const G4VProcess* aProcess1 = point1->GetProcessDefinedStep();
//    G4cout << aProcess1->GetProcessName() << G4endl;
//    G4StepPoint* point2 = step->GetPostStepPoint();
//    const G4VProcess* aProcess2 = point2->GetProcessDefinedStep();
//    G4cout << aProcess2->GetProcessName() << G4endl;
    
	// Get the interaction volume
    G4String volumeName = track -> GetVolume() -> GetName();

    targetSolid = track -> GetVolume() -> GetLogicalVolume() -> GetSolid();
    targetVisExtent = targetSolid->GetExtent();

    x_center = targetVisExtent.GetExtentCenter()[0];
    y_center = targetVisExtent.GetExtentCenter()[1];
    z_center = targetVisExtent.GetExtentCenter()[2];


    G4String volumeNum = volumeName.substr(volumeName.size()-12, 3);
    G4int volumeNumber = std::stoi(volumeNum);
    G4cout << "*** DEBUG *** volumeName: " << volumeName << " VolumeNum (str): " << volumeNum  <<  " VolumeNum (int): " << volumeNumber << G4endl;
    G4cout << "*** DEBUG *** volume position (x, y, z of the center): " << x_center << " " << y_center << " " << z_center << G4endl;
    

    G4int copyIDinX   =  volumeNumber;
    G4int copyIDinY   =  0;
    G4int icell = 0;
    
    
    
    if(CellID[copyIDinX][copyIDinY]==-1) // if no energy depositions before in this cell
	{
		// Now we create a new hit object, corresponding to the cell, and fill it with values to be stored
		DetectorHit* hit = new DetectorHit();
		hit -> SetEnergyDeposit(energyDeposit);
		hit -> SetScintID(volumeNumber);

        hit -> SetScint_x(x_center);
        hit -> SetScint_y(y_center);
        hit -> SetScint_z(z_center);

		hit -> SetTime(globaltime);
		icell = hitsCollection -> insert(hit);
        CellID[copyIDinX][copyIDinY] = icell - 1;
    }
	else // if energy depositions before in this cell
	{
		(*hitsCollection)[0] -> AddEnergyDeposit(energyDeposit);
		(*hitsCollection)[0] -> SetTime(globaltime);
	}
	return true;
}


void SensitiveDetector::EndOfEvent(G4HCofThisEvent* hitsCollEvent) 
{
//	// Method called at the end of an event.
//	if(hitsCollectionID < 0) {
//		hitsCollectionID =
//		G4SDManager::GetSDMpointer() -> GetCollectionID(hitsCollection);
//	}
//	// Finally, the hits collection just defined (with its ID and its hits) 
//	// is registered to the HCOfThisEvent of the current event. 
//	hitsCollEvent -> AddHitsCollection(hitsCollectionID, hitsCollection);
}
