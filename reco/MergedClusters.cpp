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
  }

  bool MergedClusters::DoReconstruct()
  {
    fClusters.clear(); //Clear out the old clusters from last time!

    //Make my own copy of the vector of hits as a std::list so I can remove the ones I use
    std::vector<std::pair<pers::MCCluster, std::vector<pers::MCHit>>> clusterToHit;

    const std::string fiducial = "volA3DST_PV";
    auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

    //Get a list of MCHits that are in the fiducial volume
    //Now done during hit-making
    /*std::list<pers::MCHit> hits(fHits.begin(), fHits.end());
    hits.remove_if([mat, shape](const auto& hit)
                  {
                    const auto local = ::InLocal(hit.Position.Vect(), mat);
                    double pos[] = {local.X(), local.Y(), local.Z()};
                    return shape->Contains(pos);
                  });*/

    //Tejin-like candidates (from Minerva).  
    for(auto outerHitPos = fHits.begin(); outerHitPos != fHits.end(); ++outerHitPos) //hits.begin(); outerHitPos != hits.end(); ++outerHitPos)
    {
      auto& outerHit = *outerHitPos; //Thanks for nothing TTreeReaderArray<pers::MCHit>::Iterator_t      

      //Look for an MCCluster that overlaps with this MCHit.  This nested find_if is basically a way of looping over 
      //all of the already-matched MCHits.  
      auto neighbor = std::find_if(clusterToHit.begin(), clusterToHit.end(), 
                                   [this, &outerHit](const auto& pair)
                                   {
                                     return std::find_if(pair.second.begin(), pair.second.end(), 
                                                         [this, &outerHit](const auto& innerHit)                                          
                                                         {
                                                           auto diff = outerHit.Position-innerHit.Position;
                                                           const double width = (outerHit.Width+innerHit.Width)/2.*this->fMergeDist;
                                                           return (std::fabs(diff.X()) < width 
                                                                   || std::fabs(diff.Y()) < width 
                                                                   || std::fabs(diff.Z()) < width);
                                                         }) != pair.second.end();
                                    });

      if(neighbor != clusterToHit.end()) //Merge outerHit with its' neighbor's MCCluster
      {
        auto& cluster = neighbor->first;
        cluster.Energy += outerHit.Energy;
        
        cluster.TrackIDs.insert(cluster.TrackIDs.end(), outerHit.TrackIDs.begin(), outerHit.TrackIDs.end());   
        neighbor->second.push_back(outerHit);
      }
      else //Seed a new MCCluster
      {
        pers::MCCluster seed;
        seed.Energy = outerHit.Energy;
        seed.TrackIDs = outerHit.TrackIDs;
        //seed.Position = outerHit.Position;
        clusterToHit.push_back(std::make_pair(seed, std::vector<pers::MCHit>({outerHit})));
      }
    }

    //Calculate cluster sizes
    for(auto& pair: clusterToHit)
    {
      auto& clust = pair.first;

      //Set cluster's position to the position of the first MCHit it conatins in time.  
      clust.Position = std::min_element(pair.second.begin(), pair.second.end(), [](const auto& first, const auto& second)
                                                                                { return first.Position.T() < second.Position.T(); })->Position;

      //Now that I know clust's starting position, find its size.
      const auto xExtrema = std::minmax_element(pair.second.begin(), pair.second.end(), [&clust](const auto& first, const auto& second)
                                                                                        {
                                                                                          return   (first.Position-clust.Position).X() 
                                                                                                 < (second.Position-clust.Position).X();
                                                                                        });
      clust.XWidth = std::fabs(xExtrema.second->Position.X()-xExtrema.first->Position.X())+(xExtrema.first->Width+xExtrema.second->Width)/2.;

      const auto yExtrema = std::minmax_element(pair.second.begin(), pair.second.end(), [&clust](const auto& first, const auto& second)
                                                                                        {
                                                                                          return   (first.Position-clust.Position).Y()
                                                                                                 < (second.Position-clust.Position).Y();
                                                                                        });
      clust.YWidth = std::fabs(yExtrema.second->Position.Y()-yExtrema.first->Position.Y())+(yExtrema.first->Width+yExtrema.second->Width)/2.;

      const auto zExtrema = std::minmax_element(pair.second.begin(), pair.second.end(), [&clust](const auto& first, const auto& second)
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

