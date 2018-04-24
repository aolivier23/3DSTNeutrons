//File: NeutronTOF.h
//Brief: NeutronTOF is an Analyzer that make histograms of quantities for Final State neutrons.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTreeReaderArray.h"

//EDepNeutrons includes
#include "ana/Analyzer.h"

//persistency includes
#include "persistency/MCHit.h"

//c++ includes
#include <random>

#ifndef ANA_NEUTRONTOF_H
#define ANA_NEUTRONTOF_H

namespace ana
{
  class NeutronTOF: public plgn::Analyzer
  {
    public:
      NeutronTOF(const plgn::Analyzer::Config& config);
      virtual ~NeutronTOF() = default;

    protected:
      virtual void DoAnalyze() override;

      TTreeReaderArray<pers::MCHit> fHits; //The source of MCHits to be analyzed

      TH1D* fNeutronHitTime; //Time of first hit caused by each ancestor of a FS neutron in ns
      TH2D* fNeutronTimeVersusDist; //Times of first hit caused by each ancestor of a FS neutron versus its' distance from interaction vertex.  
      TH1D* fNeutronEResidual; //((Neutron energy from TOF) - (true energy))/(true energy)
      TH1D* fNeutronTOFEnergy; //Neutron energy from TOF and distance
      //TH2D* fNeutronHitTimeVersusTrueEnergy; //Neutron hit time versus true energy
      
      TH1D* fBeta; //v/c
      TH1D* fTrueBeta; //v/c from true FS neutron initial energy
      TH1D* fBetaRes; //How well can I tell that beta is not really 1 
      TH1D* fFSNeutronEnergy; //Energies of FS neutrons "reconstructed" so that I can report an efficiency for any cuts made here

      //PRNG for smearing vertex times
      std::mt19937 fGen; //Mersenne Twister engine with period of 19937
      std::normal_distribution<double> fGaus; //Normal distribution object (Gaussian distribution)

      //Configuration parameters I want to keep around for statistics
      double fPosRes; //Position resolution for MCHits
      double fTimeRes; //Time resolution for MCHits
  };
}

#endif //ANA_NEUTRONCAND_H
