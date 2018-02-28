//File: MergedClusters.cpp
//Brief: Combines all MCHits that are adjacent to other MCHits into one big cluster.  Then, combines leftover MCHits into clusters that are 5 or fewer hit widths away.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

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

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<reco::MergedClusters>(opt::CmdLine& opts)
  {
    opts.AddKey("--hit-alg", "Name of the branch from which to read pers::MCHits.  Usually also the name of the hit-making algorithm.", "GridNeutronHits");
    opts.AddKey("--merge-dist", "Number of empty cubes over which clusters can \"jump\".  A value of 0 results in clusters of directly adjacent hits only.", "0");
  }
}

namespace reco
{
  MergedClusters::MergedClusters(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fClusters(), 
                                                                             fHits(*(config.Input), (*(config.Options))["--hit-alg"].c_str())
  {
    config.Output->Branch("MergedClusters", &fClusters);
    fMergeDist = config.Options->Get<size_t>("--merge-dist");
    fHitAlgName = (*(config.Options))["--hit-alg"];
  }

  bool MergedClusters::DoReconstruct()
  {
    fClusters.clear(); //Clear out the old clusters from last time!

    //Make my own copy of the vector of hits as a std::list so I can remove the ones I use
    std::list<pers::MCCluster> clusters;

    const std::string fiducial = "volA3DST_PV";
    auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

    //Tejin-like candidates (from Minerva).  
    for(auto outerHitPos = fHits.begin(); outerHitPos != fHits.end(); ++outerHitPos)
    {
      auto& outerHit = *outerHitPos; 

      //Prepare a new MCCluster with this hit.  I will accumulate all clusters that are within fMergeDist of this hit into seed.
      pers::MCCluster seed;
      seed.Energy = outerHit.Energy;
      seed.TrackIDs = outerHit.TrackIDs;
      seed.HitAlg = fHitAlgName;
      seed.Hits = std::vector<size_t>({outerHitPos.fIndex});

      //Look for clusters that are close to hit, meging them into seed as I go
      clusters.remove_if([&outerHit, this, &seed](const auto& cluster)
                         {
                           std::vector<pers::MCHit> hits;
                           for(const auto& pos: cluster.Hits) hits.push_back(fHits[pos]);
                           if(std::find_if(hits.begin(), hits.end(), [this, &outerHit](const auto& innerHit)
                                                                     {
                                                                       auto diff = outerHit.Position-innerHit.Position;
                                                                       const double width = (outerHit.Width+innerHit.Width)/2.*(this->fMergeDist+1.);
                                                                       return (std::fabs(diff.X()) < width
                                                                               || std::fabs(diff.Y()) < width
                                                                               || std::fabs(diff.Z()) < width);
                                                                     }) != hits.end())
                           //Merge cluster into seed and delete it
                           {
                             seed.Energy += cluster.Energy;
                             seed.TrackIDs.insert(seed.TrackIDs.end(), cluster.TrackIDs.begin(), cluster.TrackIDs.end());
                             seed.Hits.insert(seed.Hits.end(), cluster.Hits.begin(), cluster.Hits.end());
                             return true;
                           }

                           return false; 
                         });
    }

    //Calculate cluster sizes
    for(auto& clust: clusters)
    {
      //Set cluster's position to the position of the first MCHit it conatins in time.  
      //clust.Position = std::min_element(pair.second.begin(), pair.second.end(), [](const auto& first, const auto& second)
      //                                                                          { return first.Position.T() < second.Position.T(); })->Position;
      std::vector<pers::MCHit> hits;
      for(const auto& pos: clust.Hits) hits.push_back(fHits[pos]);
      
      //Set cluster's position to energy-weighted centroid
      clust.Position = std::accumulate(hits.begin(), hits.end(), TLorentzVector(0., 0., 0., 0.), 
                                       [](auto sum, const auto& hit) { return sum+hit.Position*hit.Energy; })*(1./clust.Energy);

      //Now that I know clust's starting position, find its size.
      const auto xExtrema = std::minmax_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second)
                                                                          {
                                                                            return   (first.Position-clust.Position).X() 
                                                                                   < (second.Position-clust.Position).X();
                                                                          });
      clust.XWidth = std::fabs(xExtrema.second->Position.X()-xExtrema.first->Position.X())+(xExtrema.first->Width+xExtrema.second->Width)/2.;

      const auto yExtrema = std::minmax_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second)
                                                                          {
                                                                            return   (first.Position-clust.Position).Y()
                                                                                   < (second.Position-clust.Position).Y();
                                                                          });
      clust.YWidth = std::fabs(yExtrema.second->Position.Y()-yExtrema.first->Position.Y())+(yExtrema.first->Width+yExtrema.second->Width)/2.;

      const auto zExtrema = std::minmax_element(hits.begin(), hits.end(), [&clust](const auto& first, const auto& second)
                                                                          {
                                                                            return   (first.Position-clust.Position).Z()
                                                                                   < (second.Position-clust.Position).Z();
                                                                          });
      clust.ZWidth = std::fabs(zExtrema.second->Position.Z()-zExtrema.first->Position.Z())+(zExtrema.first->Width+zExtrema.second->Width)/2.;
      fClusters.push_back(clust);
    }

    return !(fClusters.empty());
  }
  REGISTER_PLUGIN(MergedClusters, plgn::Reconstructor);
}

