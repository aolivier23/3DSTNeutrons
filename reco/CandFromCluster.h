//File: CandFromCluster.h
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates exactly one NeutronCand from each MCCluster it contains.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/NeutronCand.h"
#include "persistency/MCCluster.h"

//ROOT includes
#include "TTreeReaderArray.h"

#ifndef RECO_CANDFROMCLUSTER_H
#define RECO_CANDFROMCLUSTER_H

namespace reco
{
  class CandFromCluster: public plgn::Reconstructor
  {
    public:
      CandFromCluster(const plgn::Reconstructor::Config& config);
      virtual ~CandFromCluster() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCClusters that will be written to tree
      std::vector<pers::NeutronCand> fCands;

      //Location from which MCClusters will be read
      TTreeReaderArray<pers::MCCluster> fClusters;

      std::string fClusterAlgName; //Name of the cluster algorithm to be stitched

      //Configuration data
      double fTimeRes; //Time resolution for 3DST
      double fPosRes; //Position resolution for 3DST
  };
}

#endif //RECO_CANDFROMCLUSTER_H
