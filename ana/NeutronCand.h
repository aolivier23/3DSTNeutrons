//File: NeutronCand.h
//Brief: NeutronCand is an Analyzer that make histograms of quantities for Final State neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTreeReaderArray.h"

//EDepNeutrons includes
#include "ana/Analyzer.h"

//persistency includes
#include "persistency/MCCluster.h"

#ifndef ANA_NEUTRONCAND_H
#define ANA_NEUTRONCAND_H

namespace ana
{
  class NeutronCand: public plgn::Analyzer
  {
    public:
      NeutronCand(const plgn::Analyzer::Config& config);
      virtual ~NeutronCand() = default;

    protected:
      virtual void DoAnalyze() override;

      TTreeReaderArray<pers::MCCluster> fClusters; //The source of MCClusters to be analyzed

      const double fMinEnergy; //Energy cut used for FS neutrons

      TH1D* fCandidateEnergy; //Histogram of neutron candidate energy
      TH1D* fCandPerNeutron; //Number of candidates per FS neutron
      TH1I* fNCand; //Number of candidates per event
      TH1D* fFSNeutronEnergy; //FS neutron energies
      TH2D* fCauseEnergyVsCandEnergy; //Energy of candidates versus the FS neutrons they are descended from
      TH1D* fCandAngleWRTCause; //Angle of candidate to FS neutron's initial momentum
      TH1D* fDistFromVtx; //Distance of candidates from neutrino vertex.  Only count first candidate for each FS neutron.
  };
}

#endif //ANA_NEUTRONCAND_H
