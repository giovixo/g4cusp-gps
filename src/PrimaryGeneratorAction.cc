#include "PrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "WriteToFile.hh"

//using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    particleGun = new G4ParticleGun(1); // 1 particle per event

    // Get particle definition
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

    // Set particle properties
    particleGun->SetParticleDefinition(particle);
        // Output for debug
    std::ofstream outputFile("output.txt");

    // Read events from CSV file
    #ifdef CSV
       // rsmSource.ReadCSV("events.csv"); // no more needed as rsmSource is a global variable
    #endif
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}

void PrimaryGeneratorAction::GenerateVertex(G4double source_h, G4double source_r, G4double source_alpha, G4ThreeVector &vertex_center, G4ThreeVector& vertex) {
    G4double y = source_h + 4.*cm ; // TdDO: check the exact position of the detector surface (detector center = 10.5*cm?)
    G4double z = source_r * std::cos(source_alpha); 
    G4double x = source_r * std::sin(source_alpha);
    vertex_center.setX(x);
    vertex_center.setY(y);
    vertex_center.setZ(z);

    G4double beam_size1{10.*cm}, beam_size2{10.*cm};
    G4double x_move = -beam_size1/2. + beam_size1*G4UniformRand();
    G4double y_move = 0;
    G4double z_move = -beam_size2/2. + beam_size2*G4UniformRand();
    G4ThreeVector translation(x_move, y_move, z_move);

    vertex = vertex_center + translation;
}

void PrimaryGeneratorAction::GenerateDirection(G4double detector_h, G4ThreeVector vertex, G4ThreeVector& direction) {
    // Define the center of the detector surface
    G4ThreeVector detector_center;
    detector_center.setX(0.);
    detector_center.setY(detector_h);
    detector_center.setZ(0.);

    direction = detector_center - vertex;

    // Normalize the direction vector
    direction = direction.unit();
}

void PrimaryGeneratorAction::GenerateSquare(G4double number) {
    G4cout << "Hello " << number * number << G4endl;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {   
    
    // *** Step 0: Get the event 

    static size_t currentEventIndex = 0;

    #ifdef CSV
       const std::vector<EventData>& events = gRsmSource->GetEvents();

       if (currentEventIndex >= events.size()) {
           G4cerr << "No more events in the CSV file." << G4endl;
           return;
       }
       // This statement accesses  events[0] and then increments the index to 1, and so on.
       const EventData& event = events[currentEventIndex++];
    #endif 

    #ifndef CSV
       // This statement accesses the event 
       EventData event = gRsmSource->GenerateEvent(60.*keV, 0.01); // 0.01 photons per second
    #endif

    // *** Step 1: retrieve the event data

    G4int eventID = event.event_id; // This can be also obtained from currentEventIndex or anEvent->GetEventID()
    G4double time = event.time;
    G4double angle = event.angle;
    G4double pol_x = event.pol_x;
    G4double pol_y = event.pol_y;
    G4double pol_z = event.pol_z;
    G4double energy = event.energy;

    // *** Step 2: Generate the primary vertex using the event data

    // Set the energy of the particle
    particleGun->SetParticleEnergy(energy); 
    // Set the source vertex
    constexpr G4double source_h = 6.*cm;
    constexpr G4double source_r = 3.*cm;
    G4ThreeVector vertex_center, vertex;
    GenerateVertex(source_h, source_r, angle, vertex_center, vertex);    
    //particleGun->SetParticlePosition(vertex);
    particleGun->SetParticlePosition(vertex_center); // FOR TESTING ONLY
    // Set the source direction
    G4ThreeVector direction;
    constexpr G4double detector_h = 4.*cm; // The height of the detector surface
    GenerateDirection(detector_h, vertex_center, direction);
    particleGun->SetParticleMomentumDirection(direction);
    // Set the source polarization
    G4ThreeVector pseudo_polarization(pol_x, pol_y, pol_z);
    G4ThreeVector polarization = pseudo_polarization.cross(direction);
    particleGun->SetParticlePolarization(polarization.unit());

    particleGun->GeneratePrimaryVertex(anEvent);

    // Choose an arbitrary vector that is not parallel to the direction vector (FOR TESTING ONLY)   
    /** G4ThreeVector arbitraryVector(1.0, 0.0, 0.0);
    if (std::abs(direction.dot(arbitraryVector)) > 0.9) {
        arbitraryVector = G4ThreeVector(0.0, 1.0, 0.0);
    }
    // Use the cross product to find a perpendicular vector
    G4ThreeVector polarization = direction.cross(arbitraryVector);    
    // Normalize the polarization vector to ensure it has unit length
    polarization = polarization.unit();
    // Define the polarization vector
    particleGun->SetParticlePolarization(polarization); **/
    
    #ifdef DEBUG
    std::ostringstream strValue;
    //strValue << vertex.getX()/cm << " " << vertex.getY()/cm << " " << vertex.getZ()/cm;
    //testOutput.print(strValue.str());
    // The event data from the csv file
    //strValue << "angle: " << angle << " Energy: " << " " << energy;
    // Generate and print new events
    //EventData new_event = rsmSource.GenerateEvent();
    strValue << "Event ID: " << event.event_id << ", "
        << "Time: " << event.time << ", "
        << "Angle: " << event.angle << ", "
        << "Pol_x: " << event.pol_x << ", "
        << "Pol_y: " << event.pol_y << ", "
        << "Pol_z: " << event.pol_z << ", "
        << "Energy: " << event.energy;
    testOutput.print(strValue.str());
    #endif
}
