//File: FSNeutrons.h
//Brief: FSNeutrons is an Analyzer that make histograms of quantities for Final State neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TH1D.h"

//c++ includes

//EDepNeutrons includes
#include "ana/Analyzer.h"

#ifndef ANA_FSNEUTRONS_H
#define ANA_FSNEUTRONS_H

namespace ana
{
  class FSNeutrons: public plgn::Analyzer
  {
    public:
      FSNeutrons(const plgn::Analyzer::Config& config);
      virtual ~FSNeutrons() = default;

    protected:
      virtual void DoAnalyze() override;

      TH1D* fNeutronEnergy; //Histogram of FS neutron energy
      TH1D* fNFSNeutrons; //Number of FS neutrons above threshold in each event

      //Command line parameters
      double fEMin; //Energy required for an FS neutron to be plotted
  };
}

#endif //ANA_FSNEUTRONS_H
