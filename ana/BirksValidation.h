//File: BirksValidation.h
//Brief: BirksValidation is an Analyzer that make histograms that validate the Birks' Law implementation in edep-sim.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TH1D.h"
#include "TH2D.h"

//c++ includes

//EDepNeutrons includes
#include "ana/Analyzer.h"

#ifndef ANA_BIRKSVALIDATION_H
#define ANA_BIRKSVALIDATION_H

namespace ana
{
  class BirksValidation: public plgn::Analyzer
  {
    public:
      BirksValidation(const plgn::Analyzer::Config& config);
      virtual ~BirksValidation() = default;

    protected:
      virtual void DoAnalyze() override;

      TH2D* fVisFracVersusdEdx; //Fraction of energy that is visible versus dE/dx
      TH1D* fBirksResidual; //Difference between simulation-corrected visible energy and correction at analysis stage

      std::map<std::string, std::pair<TH2D*, TH1D*>> fPlotsPerParticle; //Validation plots for different particle types

      //Configuration data
      double fEMin; //Hit segements with energy below this value will not be plotted
  };
}

#endif //ANA_BIRKSVALIDATION_H
