//File: BirksValidation.cpp
//Brief: BirksValidation is an Analyzer that plots kinematic quantities for Final State Neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/BirksValidation.h"
#include "app/Factory.cpp"

//util includes
#include "ROOT/Base/TFileSentry.h"

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::BirksValidation>(opt::CmdLine& opts)
  {
    opts.AddKey<opt::Accumulate>("--birks-particle", "Energy deposits for particles with this name will be plotted on their own as well "
                                                     "as with all other particle types.", "");
    opts.AddKey("--birks-E-min", "Minimum energy of deposits plotted by BirksValiation Analyzer.", "0.");
  }
}*/

namespace ana
{
  BirksValidation::BirksValidation(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fEMin(config.Options["EMin"].as<double>())
  {
    fVisFracVersusdEdx = config.File->make<TH2D>("VisFracVsdEdx", "Visible Fraction of Energy versus #frac{dE}{dx};#frac{dE}{dx};Fraction Visible;"
                                                                  "Hit Segments", 1000, 0., 300, 1000, 0, 1);
    fBirksResidual = config.File->make<TH1D>("BirksResidual", "Fractional Difference Between Visible Energy and Direct Birks' Law;"
                                                              "Fractional Residual;Hit Segments", 4000, -2, 2);

    const auto partNames = config.Options["BirksParticle"].as<std::vector<std::string>>();
    for(const auto& name: partNames)
    {
      if(name == "") continue; //TODO: Fix this feature of opt::Accumulate
      fPlotsPerParticle[name] = std::make_pair(config.File->make<TH2D>((name+"VisFracVsdEdx").c_str(), ("Visible Fraction of Energy "
                                               "versus #frac{dE}{dx} for "+name+";#frac{dE}{dx};Fraction Visible;HitSegments").c_str(), 
                                               1000, 0., 300, 1000, 0, 1),
                                               config.File->make<TH1D>((name+"BirksResidual").c_str(), ("Fractional Difference Between "
                                               "Visible Energy and Direct Birks' Law for "+name+";Fractional Residual;Hit Segments").c_str(), 
                                               4000, -2, 2));
    }
  }

  void BirksValidation::DoAnalyze()
  {
    const auto& parts = fEvent->Trajectories;

    //TODO: Make these plots per detector?
    for(const auto& det: fEvent->SegmentDetectors)
    {
      for(const auto& seg: det.second)
      { 
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const float primDep = seg.GetEnergyDeposit();
        #else
        const float primDep = seg.EnergyDeposit;
        #endif
        if(primDep >= fEMin)
        {
          #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
          const double dE = seg.GetEnergyDeposit();
          const double dx = (seg.GetStop()-seg.GetStart()).Vect().Mag();
          const double parent = seg.GetPrimaryId();
          const double deposit = seg.GetSecondaryDeposit();
          #else
          const double dE = seg.EnergyDeposit;
          const double dx = (seg.Stop-seg.Start).Vect().Mag();
          const double parent = seg.PrimaryId;
          const double deposit = seg.SecondaryDeposit;
          #endif

          fVisFracVersusdEdx->Fill(dE/dx, deposit/dE);
          const double corrected = dE/(1.+dE/dx*0.126); //TODO: Set Birks' constant for other materials
          fBirksResidual->Fill((corrected-deposit)/deposit); 

          #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
          const std::string parentName(parts[parent].GetName());
          #else
          const std::string parentName = parts[parent].Name;
          #endif

          auto found = fPlotsPerParticle.find(parentName);
          if(found != fPlotsPerParticle.end())
          {
            found->second.first->Fill(dE/dx, deposit/dE);
            found->second.second->Fill((corrected-deposit)/deposit);
          }
        }
      }
    } 
  }

  REGISTER_PLUGIN(BirksValidation, plgn::Analyzer)
}
