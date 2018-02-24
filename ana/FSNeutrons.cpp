//File: FSNeutrons.cpp
//Brief: FSNeutrons is an Analyzer that plots kinematic quantities for Final State Neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/FSNeutrons.h"
#include "app/Factory.cpp"

//util includes
#include "ROOT/Base/TFileSentry.h"
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::FSNeutrons>(opt::CmdLine& opts)
  {
    opts.AddKey("--E-min", "Minimum energy for a neutron candidate to be plotted.  Should match hit-making and cluster-making algorithms", "2.0");
  }
}

namespace ana
{
  FSNeutrons::FSNeutrons(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fEMin(config.Options->Get<double>("--E-min"))
  {
    fNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of All FS Neutrons;Energy [MeV];FS Neutrons", 200, 0, 3000);
  }

  void FSNeutrons::DoAnalyze()
  {  
    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        if(part.PDGCode == 2112 && part.Momentum.E() - part.Momentum.Mag() > fEMin) fNeutronEnergy->Fill(part.Momentum.E() - part.Momentum.Mag());
      }
    }
  }

  REGISTER_PLUGIN(FSNeutrons, plgn::Analyzer);
}
