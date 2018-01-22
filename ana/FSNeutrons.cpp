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
    fNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "Energy Spectrum of Final State Neutrons;Energy [MeV];Events", 
                                             300, 0, 2000);
  }

  void FSNeutrons::DoAnalyze()
  {
    for(const auto& traj: fEvent->Trajectories)
    {
      if(traj.PDGCode == 2112) fNeutronEnergy->Fill(traj.InitialMomentum.E());
    }
  }

  REGISTER_PLUGIN(FSNeutrons, plgn::Analyzer);
}
