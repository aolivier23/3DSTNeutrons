//File: FSNeutrons.cpp
//Brief: FSNeutrons is an Analyzer that plots kinematic quantities for Final State Neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/FSNeutrons.h"
#include "app/Factory.cpp"

//util includes
#include "ROOT/Base/TFileSentry.h"

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::FSNeutrons>(opt::CmdLine& opts)
  {
    opts.AddKey("--E-min", "Minimum energy for a neutron candidate to be plotted.  Should match hit-making and cluster-making algorithms", "1.5");
  }
}*/

namespace ana
{
  FSNeutrons::FSNeutrons(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fEMin(config.Options["EMin"].as<double>())
  {
    fNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of All FS Neutrons;Energy [MeV];FS Neutrons", 200, 0, 3000);
    fNFSNeutrons = config.File->make<TH1D>("NFSNeutrons", ("Number of FS Neutrons Above "+std::to_string(fEMin)+" MeV;FS Neutrons;Events").c_str(), 
                                           10, 0, 10);
  }

  void FSNeutrons::DoAnalyze()
  {  
    size_t nFSNeutrons = 0;
    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const double KE = part.GetMomentum().E() - part.GetMomentum().Mag();
        const int pdg = part.GetPDGCode();
        #else
        const double KE = part.Momentum.E() - part.Momentum.Mag();
        const int pdg = part.PDGCode;
        #endif

        if(pdg == 2112 && KE > fEMin) 
        {
          fNeutronEnergy->Fill(KE);
          ++nFSNeutrons;
        }
      }
    }
    fNFSNeutrons->Fill(nFSNeutrons);
  }

  REGISTER_PLUGIN(FSNeutrons, plgn::Analyzer);
}
