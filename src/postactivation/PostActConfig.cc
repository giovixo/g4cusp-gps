#include "PostActConfig.hh"

#include "G4GenericMessenger.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4NuclideTable.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumeStore.hh"

#include <regex>

G4String   PostActConfig::fIsotopeName = "";
PostActIon PostActConfig::fIon;
G4String   PostActConfig::fG4IonName = "";
G4String   PostActConfig::fVolumeName = "";
G4String   PostActConfig::fOutputPrefix = "postact";
G4double   PostActConfig::fMinHalfLife = 1*us;
G4double   PostActConfig::fTimeWindow = 10*us;


G4GenericMessenger* PostActConfig::CreateMessenger()
{
    // The messenger owns the object whose methods it calls; the values are static
    static PostActConfig config;
    auto messenger = new G4GenericMessenger(&config, "/postact/", "Post-activation decays");

    messenger->DeclareMethod("isotope", &PostActConfig::SetIsotope,
                             "Nuclide to decay, e.g. Na22, Co60[58.590], Ta178[0.000X]")
        .SetParameterName("name", false)
        .SetStates(G4State_Idle)          // the ion is created here: needs an initialised kernel
        .SetToBeBroadcasted(false);
    messenger->DeclareMethod("volume", &PostActConfig::SetVolume,
                             "Physical volume where the nuclides are placed")
        .SetParameterName("name", false)
        .SetToBeBroadcasted(false);
    messenger->DeclareProperty("output", fOutputPrefix, "Prefix of the output files")
        .SetParameterName("prefix", false)
        .SetToBeBroadcasted(false);
    messenger->DeclarePropertyWithUnit("minHalfLife", "s", fMinHalfLife,
                                       "Daughters with T1/2 >= this value are killed before they decay")
        .SetParameterName("t", false)
        .SetRange("t>0")
        .SetToBeBroadcasted(false);
    messenger->DeclarePropertyWithUnit("timeWindow", "s", fTimeWindow,
                                       "Deposits later than this after the decay are dropped")
        .SetParameterName("t", false)
        .SetRange("t>0")
        .SetToBeBroadcasted(false);
    return messenger;
}


G4bool PostActConfig::ParseIonName(const G4String& name, PostActIon& ion)
{
    static const std::regex re(R"(^([A-Z][a-z]?)(\d+)(?:\[(\d+(?:\.\d*)?)([XYZUVW])?\])?$)");
    std::smatch m;
    if (!std::regex_match(name, m, re)) return false;

    const G4int Z = G4lrint(G4NistManager::Instance()->GetZ(m[1].str()));
    if (Z <= 0) return false;

    ion.Z = Z;
    ion.A = std::stoi(m[2].str());
    ion.E = m[3].matched ? std::stod(m[3].str())*keV : 0.;
    ion.flb = m[4].matched ? m[4].str()[0] : '\0';
    return ion.A >= Z;
}


void PostActConfig::SetIsotope(const G4String& name)
{
    PostActIon ion;
    if (!ParseIonName(name, ion))
    {
        G4ExceptionDescription msg;
        msg << "Invalid nuclide name \"" << name << "\" (expected e.g. Na22, Co60[58.590], Ta178[0.000X])";
        G4Exception("PostActConfig::SetIsotope", "PostAct001", FatalErrorInArgument, msg);
        return;
    }

    const G4Ions::G4FloatLevelBase flb = (ion.flb == '\0') ? G4Ions::G4FloatLevelBase::no_Float
                                                          : G4Ions::FloatLevelBase(ion.flb);

    // The pipeline names give E rounded to 1 eV (e.g. I130[39.953] for the 39.9525 keV level),
    // while the nuclide table matches within +-0.5 eV: look for the level within the rounding
    // and use its exact energy
    if (ion.E > 0.)
    {
        const G4IsotopeProperty* level = nullptr;
        for (const G4double dE : {0., -0.4999*eV, 0.4999*eV})
        {
            level = G4NuclideTable::GetNuclideTable()->GetIsotope(ion.Z, ion.A, ion.E + dE, flb);
            if (level != nullptr) break;
        }
        if (level == nullptr)
        {
            G4ExceptionDescription msg;
            msg << "Isotope " << name << ": no level at " << ion.E/keV << " keV in the Geant4 nuclide "
                << "table (G4IonTable would create an artificial one)";
            G4Exception("PostActConfig::SetIsotope", "PostAct019", FatalException, msg);
            return;
        }
        ion.E = level->GetEnergy();
    }

    // Create the ion on the master, so that all threads decay the same Geant4 level
    const auto* g4ion = dynamic_cast<const G4Ions*>(G4IonTable::GetIonTable()->GetIon(ion.Z, ion.A, ion.E, flb));

    // Never decay a nuclide different from the requested one
    G4ExceptionDescription msg;
    msg << "Isotope " << name << ": ";
    if (g4ion == nullptr)
    {
        msg << "G4IonTable::GetIon returned no ion";
        G4Exception("PostActConfig::SetIsotope", "PostAct012", FatalException, msg);
        return;
    }
    if (std::abs(g4ion->GetExcitationEnergy() - ion.E) > 1*eV)
    {
        msg << "requested excitation energy " << ion.E/keV << " keV, but Geant4 gave "
            << g4ion->GetParticleName();
        G4Exception("PostActConfig::SetIsotope", "PostAct013", FatalException, msg);
        return;
    }
    if (ion.E > 0.)
    {
        if (g4ion->GetFloatLevelBase() != flb)
        {
            msg << "requested floating level '" << ion.flb << "', but Geant4 gave " << g4ion->GetParticleName();
            G4Exception("PostActConfig::SetIsotope", "PostAct013", FatalException, msg);
            return;
        }
    }
    else if (g4ion->GetFloatLevelBase() != flb)
    {
        // E = 0: Geant4 has a single level per nuclide (see G4IonTable::FindIon)
        msg << "Geant4 has only one level at E = 0 for this nuclide; it decays as "
            << g4ion->GetParticleName() << ", as a daughter " << name << " would in Geant4";
        G4Exception("PostActConfig::SetIsotope", "PostAct022", JustWarning, msg);
    }

    fIsotopeName = name;
    fIon = ion;
    fG4IonName = g4ion->GetParticleName();
}


void PostActConfig::SetVolume(const G4String& name)
{
    if (G4PhysicalVolumeStore::GetInstance()->GetVolume(name, false) == nullptr)
    {
        G4ExceptionDescription msg;
        msg << "Physical volume \"" << name << "\" not found";
        G4Exception("PostActConfig::SetVolume", "PostAct002", FatalErrorInArgument, msg);
        return;
    }
    fVolumeName = name;
}
