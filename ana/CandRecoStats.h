//File: CandRecoStats.h
//Brief: CandRecoStats is an Analyzer that make histograms of quantities for Final State neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTreeReaderArray.h"

//EDepNeutrons includes
#include "ana/Analyzer.h"

//persistency includes
#include "persistency/NeutronCand.h"

#ifndef ANA_CANDRECOSTATS_H
#define ANA_CANDRECOSTATS_H

namespace ana
{
  class CandRecoStats: public plgn::Analyzer
  {
    public:
      CandRecoStats(const plgn::Analyzer::Config& config);
      virtual ~CandRecoStats() = default;

    protected:
      virtual void DoAnalyze() override;

      TTreeReaderArray<pers::NeutronCand> fCands; //The source of NeutronCands to analyze

      const double fMinEnergy; //Energy cut used for FS neutrons

      TH1D* fCandidateEnergy; //Histogram of neutron candidate energy
      TH1D* fCandPerNeutron; //Number of candidates per FS neutron
      TH1I* fNCand; //Number of candidates per event
      TH1D* fFSNeutronEnergy; //FS neutron energies
      TH2D* fCauseEnergyVsCandEnergy; //Energy of candidates versus the FS neutrons they are descended from
      TH1D* fCandAngleWRTCause; //Angle of candidate to FS neutron's initial momentum
      TH1D* fDistFromVtx; //Distance of candidates from neutrino vertex.  Only count first candidate for each FS neutron.
      TH2D* fCandPerNeutronVsNeutronKE; //Number of candidates versus neutron KE.  Filled once per FS neutron with candidates.
      TH2D* fAngleVsDistFromVtx; //Angle of candidate w.r.t. FS neutron versus distance from vertex to show why angle distribution has 
                                 //peak at 0
      TH1D* fNeutronsPerCand; //Number of neutrons assigned to each candidate
      //TH1D* fLostNeutronE; //When multiple neutrons are grouped into one candidate, what are the energies of the neutrons that 
                           //are not the most energetic?
  };
}

#endif //ANA_CANDRECOSTATS_H
