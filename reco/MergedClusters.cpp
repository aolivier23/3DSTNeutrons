//File: MergedClusters.cpp
//Brief: Combines all MCHits that are adjacent to other MCHits into one big cluster.  Then, combines leftover MCHits into clusters that are 5 or fewer hit widths away.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/MergedClusters.h"
#include "reco/alg/GeoFunc.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

//c++ includes
#include <numeric> //std::accumulate got moved here in c++14

namespace reco
{
  MergedClusters::MergedClusters(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fClusters(), 
                                                                             fHits(*(config.Input), 
                                                                                   config.Options["HitAlg"].as<std::string>().c_str())
  {
    config.Output->Branch("MergedClusters", &fClusters);
    fMergeDist = config.Options["MergeDist"].as<size_t>();
    fHitAlgName = config.Options["HitAlg"].as<std::string>();
  }

  bool MergedClusters::DoReconstruct()
  {
    fClusters.clear(); //Clear out the old clusters from last time!

    //Make my own copy of the vector of hits as a std::list so I can remove the ones I use
    std::list<std::pair<pers::MCCluster, std::vector<pers::MCHit>>> clusterToHits;

    //Tejin-like candidates (from Minerva).  
    for(auto outerHitPos = fHits.begin(); outerHitPos != fHits.end(); ++outerHitPos)
    {
      auto& outerHit = *outerHitPos; 

      //Prepare a new MCCluster with this hit.  I will accumulate all clusters that are within fMergeDist of this hit into seed.
      std::pair<pers::MCCluster, std::vector<pers::MCHit>> seed;
      seed.first.Energy = outerHit.Energy;
      seed.first.TrackIDs = outerHit.TrackIDs;
      seed.second.push_back(outerHit);

      //Look for clusters that are close to hit, meging them into seed as I go
      clusterToHits.remove_if([&outerHit, this, &seed](const auto& pair)
                              {
                                const auto& cluster = pair.first;
                                const auto& hits = pair.second;
                                if(std::find_if(hits.begin(), hits.end(), [this, &outerHit](const auto& innerHit)
                                                                          {
                                                                            auto diff = outerHit.Position-innerHit.Position;
                                                                            const double width = (outerHit.Width+innerHit.Width)/2.*(this->fMergeDist+1.001);
                                                                            return (std::fabs(diff.X()) < width
                                                                                    && std::fabs(diff.Y()) < width
                                                                                    && std::fabs(diff.Z()) < width);
                                                                          }) != hits.end())
                                //Merge cluster into seed and delete it
                                {
                                  seed.first.Energy += cluster.Energy;
                                  seed.first.TrackIDs.insert(seed.first.TrackIDs.end(), cluster.TrackIDs.begin(), cluster.TrackIDs.end());
                                  seed.second.insert(seed.second.end(), pair.second.begin(), pair.second.end());
                                  return true;
                                }
     
                                return false; 
                              });
    
      clusterToHits.push_back(seed); //Put this cluster into the list of all clusters
    }

    //Get vertex position for deciding which MCHit is the closest to vertex.
    #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
    const auto& vertPos = fEvent->Primaries.front().GetPosition(); //TODO: What should I do if there are multiple vertices?
    #else
	const auto& vertPos = fEvent->Primaries.front().Position;
    #endif

    //Calculate cluster sizes
    for(auto& pair: clusterToHits)
    {
      //Set cluster's position to the position of the first MCHit it conatins in time.  
      //clust.Position = std::min_element(pair.second.begin(), pair.second.end(), [](const auto& first, const auto& second)
      //                                                                          { return first.Position.T() < second.Position.T(); })->Position;
      auto& clust = pair.first;
      const auto& hits = pair.second;
      
      //Set cluster's position to energy-weighted centroid
      clust.Position = std::accumulate(hits.begin(), hits.end(), TLorentzVector(0., 0., 0., 0.), 
                                       [](auto sum, const auto& hit) { return sum+hit.Position*hit.Energy; })*(1./clust.Energy);

      //Set cluster's starting position to hit closest to the vertex.  
      clust.FirstPosition = std::min_element(hits.begin(), hits.end(), [&vertPos](const auto& first, const auto& second)
                                                                       {
                                                                         return (first.Position-vertPos).Vect().Mag() < 
                                                                                (second.Position-vertPos).Vect().Mag();
                                                                       })->Position;

      //Now that I know clust's starting position, find its size.
      const auto xMax = std::max_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second) 
                                                                   { return std::fabs(first.Position.X() - clust.Position.X())
                                                                          < std::fabs(second.Position.X() - clust.Position.X()); });
      clust.XWidth = (xMax != hits.end())?2.*std::fabs(xMax->Position.X() - clust.Position.X())+xMax->Width:-1.;

      const auto yMax = std::max_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second)
                                                                   {
                                                                     return std::fabs(first.Position.Y() - clust.Position.Y())                                                                                           < std::fabs(second.Position.Y() - clust.Position.Y());
                                                                   });
      clust.YWidth = (yMax != hits.end())?2.*std::fabs(yMax->Position.Y() - clust.Position.Y())+yMax->Width:-1.;

      const auto zMax = std::max_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second)
                                                                   {
                                                                     return std::fabs(first.Position.Z() - clust.Position.Z())    
                                                                          < std::fabs(second.Position.Z() - clust.Position.Z());
                                                                   });
      clust.ZWidth = (zMax != hits.end())?2.*std::fabs(zMax->Position.Z() - clust.Position.Z())+zMax->Width:-1.;
      fClusters.push_back(clust);
    }

    return !(fClusters.empty());
  }
  REGISTER_PLUGIN(MergedClusters, plgn::Reconstructor)
}

