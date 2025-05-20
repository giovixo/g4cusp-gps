#include "RsmSource.hh"
#include <fstream>
#include <sstream>
#include <iostream>

void RsmSource::ReadCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    // Clear any existing data
    events.clear();

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        EventData event;
        std::string value;

        std::getline(ss, value, ',');
        event.event_id = std::stoi(value);
        std::getline(ss, value, ',');
        event.time = std::stod(value);
        std::getline(ss, value, ',');
        event.angle = std::stod(value);
        std::getline(ss, value, ',');
        event.pol_x = std::stod(value);
        std::getline(ss, value, ',');
        event.pol_y = std::stod(value);
        std::getline(ss, value, ',');
        event.pol_z = std::stod(value);
        std::getline(ss, value, ',');
        event.energy = std::stod(value);

        events.push_back(event);
    }

    if (file.is_open()) {
        file.close();
    }

    std::cout << "Read " << events.size() << " events from " << filename << std::endl;
}

void RsmSource::PrintCSV() const {
    std::cout << "Event ID, Time, Angle, Pol_x, Pol_y, Pol_z, Energy" << std::endl;
    for (const auto& event : events) {
        std::cout << event.event_id << ", "
                  << event.time << ", "
                  << event.angle << ", "
                  << event.pol_x << ", "
                  << event.pol_y << ", "
                  << event.pol_z << ", "
                  << event.energy << std::endl;
    }
}

EventData RsmSource::GenerateEvent(G4double energy_input, G4double photon_rate) { 
    EventData event;
    // Assign a unique event ID
    event.event_id = next_event_id++;

    // Angular frequency of the rotating source
    // Assuming a rotation period of 100 seconds
    const double omega = 0.01 * (2.0 * M_PI);

    // Generate a random time interval
    double time_interval = G4RandExponential::shoot(1.0 / photon_rate);
    event.time += time_interval;
    event.angle = omega * event.time;
    event.pol_x = std::cos(event.angle);
    event.pol_y = 0.0;
    event.pol_z = std::sin(event.angle);
    event.energy = energy_input;

    return event;
}