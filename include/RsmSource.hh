#ifndef RSMSOURCE_HH
#define RSMSOURCE_HH

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include "Randomize.hh" // GEANT4 random number generator
#include "globals.hh" // GEANT4 global definitions (like G4double and G4int)
#include "G4SystemOfUnits.hh" // GEANT4 units 

struct EventData {
    int event_id;
    double time;
    double angle;
    double pol_x;
    double pol_y;
    double pol_z;
    double energy;
};

class RsmSource {
public:
    RsmSource() = default;
    ~RsmSource() = default;

    void ReadCSV(const std::string& filename);
    const std::vector<EventData>& GetEvents() const { return events; }
    EventData GenerateEvent(G4double energy_input, G4double photon_rate); // Method to generate an EventData
    void PrintCSV() const; // Method to print the CSV data
    G4int GetEventID() const { return next_event_id-1; } // Method to get the next event ID
    G4double GetTime () const { return events[next_event_id-1].time; } // Method to get the time of the next event

private:
    std::vector<EventData> events;
    int next_event_id = 0; // Keep track of the next event ID
};

#endif // RSMSOURCE_HH