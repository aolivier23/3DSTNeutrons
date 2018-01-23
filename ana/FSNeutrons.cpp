//File: FSNeutrons.cpp
//Brief: FSNeutrons is an Analyzer that plots kinematic quantities for Final State Neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/FSNeutrons.h"
#include "app/Factory.cpp"

//util includes
#include "ROOT/Base/TFileSentry.h"

namespace ana
{
  FSNeutrons::FSNeutrons(const plgn::Analyzer::Config& config): plgn::Analyzer(config)
  {
    fNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "Energies of All FS Neutrons;Energy [MeV];FS Neutrons", 300, 0, 2000);
  }

  void FSNeutrons::DoAnalyze()
  {  
    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        if(part.PDGCode == 2112) fNeutronEnergy->Fill(part.Momentum.E());
      }
    }
  }

  REGISTER_PLUGIN(FSNeutrons, plgn::Analyzer);
}
