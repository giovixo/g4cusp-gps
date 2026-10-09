#include "PostActConfig.hh"

#include "G4GenericMessenger.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumeStore.hh"

#include <regex>

G4String   PostActConfig::fIsotopeName = "";
PostActIon PostActConfig::fIon;
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
    fIsotopeName = name;
    fIon = ion;
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
