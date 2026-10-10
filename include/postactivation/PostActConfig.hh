#ifndef PostActConfig_h
#define PostActConfig_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4GenericMessenger;

// Settings of cusp-postactivation, shared by all threads.
// They are set on the master by the /postact/ commands, before /run/beamOn, and are
// read-only during a run. One run = one (volume, isotope) pair.
//
//   /postact/isotope   Ta178[0.000X]     nuclide to decay (step-3 naming: Sym A [E keV + floating letter])
//                                        The ion is created on the master when the command is given.
//                                        At E = 0, Geant4 creates only one level per nuclide (the first
//                                        one requested, or its preferred floating level): if it differs
//                                        from the requested one, the nuclide decays as that level, as a
//                                        daughter would in Geant4, and a warning is printed.
//   /postact/volume    PV-Absorber_006   physical volume where the nuclides are placed
//   /postact/output    postact           prefix of the output files
//   /postact/minHalfLife 1 us            daughters with T1/2 >= this are killed before decaying
//   /postact/timeWindow 10 us            deposits later than this after the decay are dropped

struct PostActIon
{
    G4int    Z = 0;
    G4int    A = 0;
    G4double E = 0.;        // excitation energy (Geant4 internal units)
    char     flb = '\0';    // floating-level letter (X, Y, Z, U, V, W) or '\0'
};

class PostActConfig
{
public:
    // Creates the /postact/ commands. Call it once, on the master thread (see main).
    static G4GenericMessenger* CreateMessenger();

    // Parses a step-3 nuclide name, e.g. "Na22", "Co60[58.590]", "Ta178[0.000X]".
    // Returns false if the name is not valid.
    static G4bool ParseIonName(const G4String& name, PostActIon& ion);

    static const G4String& IsotopeName()  { return fIsotopeName; }
    static const PostActIon& Ion()        { return fIon; }
    static const G4String& G4IonName()    { return fG4IonName; }   // name of the Geant4 ion that decays
    static const G4String& VolumeName()   { return fVolumeName; }
    static const G4String& OutputPrefix() { return fOutputPrefix; }
    static G4double MinHalfLife()         { return fMinHalfLife; }
    static G4double TimeWindow()          { return fTimeWindow; }

private:
    void SetIsotope(const G4String& name);
    void SetVolume(const G4String& name);

    static G4String   fIsotopeName;
    static PostActIon fIon;
    static G4String   fG4IonName;
    static G4String   fVolumeName;
    static G4String   fOutputPrefix;
    static G4double   fMinHalfLife;
    static G4double   fTimeWindow;
};

#endif
