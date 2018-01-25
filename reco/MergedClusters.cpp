//File: MergedClusters.cpp
//Brief: Combines all MCHits that are adjacent to other MCHits into one big cluster.  Then, combines leftover MCHits into clusters that are 5 or fewer hit widths away.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/MergedClusters.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

namespace
{
  //TODO: If I am going to reuse these, put them into a class or at least a common header
  //Return the product of the matrices from the node with a volume called name with all of its' ancestors.
  TGeoMatrix* findMat(const std::string& name, TGeoNode& parent)
  {
    if(std::string(parent.GetVolume()->GetName()) == name) return parent.GetMatrix();
    auto children = parent.GetNodes();
    for(auto child: *children)
    {
      auto node = (TGeoNode*)child;
      auto result = findMat(name, *node);
      if(result)
      {
        auto retVal = new TGeoHMatrix(*(parent.GetMatrix()));
        retVal->Multiply(result);
        return retVal;
      }
    }
    return nullptr;
  }

  //Return a TVector3 in a different coordinate system
  TVector3 InLocal(const TVector3& pos, TGeoMatrix* mat)
  {
    double master[3] = {}, local[3] = {};
    pos.GetXYZ(master);
    mat->MasterToLocal(master, local);
    return TVector3(local[0], local[1], local[2]);
  }
}

namespace reco
{
  MergedClusters::MergedClusters(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fClusters(), fHits(*(config.Input), "NoGridNeutronHits")
  {
    config.Output->Branch("MergedClusters", &fClusters);
  }

  bool MergedClusters::DoReconstruct()
  {
    fClusters.clear(); //Clear out the old clusters from last time!

    //Make my own copy of the vector of hits as a std::list so I can remove the ones I use
    std::vector<std::pair<pers::MCCluster, std::vector<pers::MCHit>>> clusterToHit;

    const std::string fiducial = "volA3DST_PV";
    auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

    //Get a list of MCHits that are in the fiducial volume
    std::list<pers::MCHit> hits(fHits.begin(), fHits.end());
    hits.remove_if([mat, shape](const auto& hit)
                  {
                    const auto local = ::InLocal(hit.Position.Vect(), mat);
                    double pos[] = {local.X(), local.Y(), local.Z()};
                    return shape->Contains(pos);
                  });

    //Tejin-like candidates (from Minerva).  
    for(auto outerHitPos = hits.begin(); outerHitPos != hits.end(); ++outerHitPos)
    {
      auto& outerHit = *outerHitPos; //Thanks for nothing TTreeReaderArray<pers::MCHit>::Iterator_t      

      //Look for an MCCluster that overlaps with this MCHit.  This nested find_if is basically a way of looping over 
      //all of the already-matched MCHits.  
      auto neighbor = std::find_if(clusterToHit.begin(), clusterToHit.end(), 
                                   [&outerHit](const auto& pair)
                                   {
                                     return std::find_if(pair.second.begin(), pair.second.end(), 
                                                         [&outerHit](const auto& innerHit)                                          
                                                         {
                                                           auto diff = outerHit.Position-innerHit.Position;
                                                           const double width = outerHit.Width+innerHit.Width;
                                                           return (std::fabs(diff.X()) < width/2. 
                                                                   || std::fabs(diff.Y()) < width/2. 
                                                                   || std::fabs(diff.Z()) < width/2.);
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

