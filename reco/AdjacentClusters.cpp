//File: AdjacentClusters.cpp
//Brief: Finds MCHits that are indirectly adjacent to each other and combines them into clusters.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/AdjacentClusters.h"

namespace reco
{
  AdjacentClusters::AdjacentClusters(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fClusters(), fHits(*(config.Input), "NeutronHits")
  {
    config.Output->Branch("AdjacentClusters", &fClusters);
  }

  bool AdjacentClusters::DoReconstruct()
  {
    fClusters.clear(); //Remove clusters from previous events!

    //Make my own copy of the vector of hits as a std::list so I can remove the ones I use
    std::list<pers::MCHit> hits(fHits.begin(), fHits.end());

    while(hits.size() > 0)
    {
      const auto& seed = *(hits.begin());
      pers::MCCluster clust;
      clust.Position = seed.Position;
      clust.XWidth = seed.Width;
      clust.YWidth = seed.Width;
      clust.ZWidth = seed.Width;

      //Make MCClusters from all MCHits that are indirectly adjacent to seed
      hits.remove_if([&seed, &clust](const auto& hit)
                     {
                       const auto diff = clust.Position-hit.Position;
                       const double adjacent = 5.+0.01;
                       if(std::fabs(diff.X()) < adjacent*seed.Width && std::fabs(diff.Y()) < adjacent*seed.Width && std::fabs(diff.Z()) < adjacent*seed.Width)
                       {
                         clust.Energy += hit.Energy;
                         clust.TrackIDs.insert(clust.TrackIDs.end(), hit.TrackIDs.begin(), hit.TrackIDs.end());
                         
                         //Update cluster size
                         if(std::fabs(diff.X()) > clust.XWidth) clust.XWidth = std::fabs(diff.X());
                         if(std::fabs(diff.Y()) > clust.YWidth) clust.YWidth = std::fabs(diff.Y());
                         if(std::fabs(diff.Z()) > clust.ZWidth) clust.ZWidth = std::fabs(diff.Z()); 

                         return true;
                       }
                       return false;
                     });
      fClusters.push_back(clust);
    }

    return !(fClusters.empty());
  }
  REGISTER_PLUGIN(AdjacentClusters, plgn::Reconstructor)
}

