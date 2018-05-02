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
      TH2D* fCandPerNeutronVsNeutronKE; //Number of candidates versus neutron KE.  Filled once per FS neutron with candidates.
      TH2D* fAngleVsDistFromVtx; //Angle of candidate w.r.t. FS neutron versus distance from vertex to show why angle distribution has 
                                 //peak at 0
      TH1D* fNeutronsPerCand; //Number of neutrons assigned to each candidate

      //Plots to investigate a likelihood-based fitting procedure.  
      TH2D* fClusterNumVsEDep; //Cluster number in time-ordering versus energy deposited in cluster
      TTree* fLikelihoodTree; //TTree with one entry for each cluster that contains information about how that cluster relates to other 
                              //clusters for its' parent true neutron
                              
      //Variables for fLikelihoodTree
      unsigned int fClusterNumber; //Position in time-ordered list of clusters for parent
      unsigned int fClustersFromEnd; //Distance from end of time-ordered list of clusters for parent
      float fDeltaAngle; //Bend in trajectory over this cluster
      float fEDep; //Energy deposited in this cluster
      float fELeft; //Energy left to be deposited by this neutron
      float fEFromTOF; //Energy from TOF for the first cluster of this cluster's parent neutron
      float fDistFromPrev; //Distance between this cluster and the next cluster
      float fDeltaT; //Time difference between clusters
      float fTrueE; //True energy of the neutron that produced a cluster
  };
}

#endif //ANA_NEUTRONCAND_H
