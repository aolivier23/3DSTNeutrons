//File: AdjacentClusters.h
//Brief: A Reconstructor that reads in a TG4Event from edepsim and creates MCHits from energy deposits that were 
//       created by ancestors of FS neutrons. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/MCHit.h"
#include "persistency/MCCluster.h"

//ROOT includes
#include "TTreeReaderArray.h"

#ifndef RECO_ADJACENTCLUSTERS_H
#define RECO_ADJACENTCLUSTERS_H

namespace reco
{
  class AdjacentClusters: public plgn::Reconstructor
  {
    public:
      AdjacentClusters(const plgn::Reconstructor::Config& config);
      virtual ~AdjacentClusters() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCClusters that will be written to tree
      std::vector<pers::MCCluster> fClusters;

      //Location from which MCHits will be read
      TTreeReaderArray<pers::MCHit> fHits;
  };
}

#endif //RECO_ADJACENTCLUSTERS_H
