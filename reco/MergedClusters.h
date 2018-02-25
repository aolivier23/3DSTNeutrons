//File: MergedClusters.h
//Brief: A Reconstructor that reads in a TG4Event from edepsim and creates MCHits from energy deposits that were 
//       created by ancestors of FS neutrons. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/MCHit.h"
#include "persistency/MCCluster.h"

//ROOT includes
#include "TTreeReaderArray.h"

#ifndef RECO_MERGEDCLUSTERS_H
#define RECO_MERGEDCLUSTERS_H

namespace reco
{
  class MergedClusters: public plgn::Reconstructor
  {
    public:
      MergedClusters(const plgn::Reconstructor::Config& config);
      virtual ~MergedClusters() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCClusters that will be written to tree
      std::vector<pers::MCCluster> fClusters;

      //Location from which MCHits will be read
      TTreeReaderArray<pers::MCHit> fHits;

      size_t fMergeDist; //Number of empty cubes over which clusters can "jump".  A value of 0 means cubes must be adjacent to form 
                         //clusters.  
  };
}

#endif //RECO_MERGEDCLUSTERS_H
