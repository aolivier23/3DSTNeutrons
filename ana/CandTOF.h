//File: CandTOF.h
//Brief: CandTOF is an Analyzer that measures the performance of NeutronCands in terms of how much neutron energy they reconstruct from TOF.  
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
#include "persistency/MCCluster.h"

//c++ includes
#include <random>

#ifndef ANA_CANDTOF_H
#define ANA_CANDTOF_H

namespace ana
{
  class CandTOF: public plgn::Analyzer
  {
    public:
      CandTOF(const plgn::Analyzer::Config& config);
      virtual ~CandTOF() = default;

    protected:
      virtual void DoAnalyze() override;

      TTreeReaderArray<pers::NeutronCand> fCands; //The source of NeutronCands to be analyzed
      TTreeReaderArray<pers::MCCluster> fClusters; //The source of MCClusters referred to by NeutronCands

      TH1D* fNeutronHitTime; //Time of first hit caused by each ancestor of a FS neutron in ns
      TH2D* fNeutronTimeVersusDist; //Times of first hit caused by each ancestor of a FS neutron versus its' distance from interaction vertex.  
      TH1D* fNeutronEResidual; //((Neutron energy from TOF) - (true energy))/(true energy)
      TH1D* fCandTOFEnergy; //Neutron energy from TOF and distance
      
      TH1D* fBeta; //v/c
      TH1D* fTrueBeta; //v/c from true FS neutron initial energy
      TH1D* fBetaRes; //How well can I tell that beta is not really 1 
      TH1D* fFSNeutronEnergy; //Energies of FS neutrons "reconstructed" so that I can report an efficiency for any cuts made here

      TH1D* fTotalEResidual; //((total energy in event) - (total energy from TOF))/(total energy in event).  How well does TOF "reconstruct" 
                             //neutron energy?
      //TH1D* fTOFELost; //Energy from TOF lost due to combining FS neutrons into the same candidate

      //PRNG for smearing vertex times
      std::mt19937 fGen; //Mersenne Twister engine with period of 19937
      std::normal_distribution<double> fGaus; //Normal distribution object (Gaussian distribution)

      //Configuration parameters I want to keep around for statistics
      double fPosRes; //Position resolution for MCHits
      double fTimeRes; //Time resolution for MCHits
  };
}

#endif //ANA_CANDTOF_H
