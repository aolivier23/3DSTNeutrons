//File: GridNeutronHits.cpp
//Brief: A Reconstructor that makes MCHits from energy deposits produced by ancestors of FS neutrons.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Base/exception.h"
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"
#include "IO/Option/runtime/Exists.h"

//edepsim includes
#include "TG4Trajectory.h"
#include "TG4HitSegment.h"

//ROOT includes
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"

//EdepNeutrons includes
#include "reco/GridNeutronHits.h"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"
#include "reco/alg/GeoFunc.h"
#include "alg/TruthFunc.h"

//c++ includes
#include <set>

namespace
{
  void RecursiveRemove(const reco::GridHits::Triple& pos, std::map<reco::GridHits::Triple, std::list<reco::GridHits::Triple>>& passed, 
                       std::map<reco::GridHits::Triple, std::list<reco::GridHits::Triple>>& failed)
  {
    //TODO: Failure in comparison operator -> I'm corrupting the map's memory somewhere?  Recursion really is a bad strategy here to begin with.  Look for something else.  
    auto found = failed.find(pos);
    if(found == failed.end()) found = passed.find(pos);
    else failed.erase(pos);
    if(found == passed.end()) return; //Nothing to do?
    else passed.erase(pos); 
                                                                                                                     
    for(const auto& neighbor: found->second) 
    {
      RecursiveRemove(neighbor, passed, failed); 
    }
  }
}

namespace reco
{
  GridNeutronHits::GridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), 
                                                                               fHitAlg(config.Options["CubeSize"].as<double>(), 
                                                                                       config.Options["AfterBirks"].as<bool>(), 
                                                                                       config.Options["TimeRes"].as<double>())
  {
    config.Output->Branch("GridNeutronHits", &fHits);
    
    fEMin = config.Options["EMin"].as<double>();
    fNeighborDist = config.Options["NeighborCut"].as<size_t>();
  } 

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool GridNeutronHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Get geometry information about this detector
    const std::string fiducial = "volA3DST_PV";
    auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();
                                                                                                                         
    //Form MCHits from all remaining hit segments
    //First, create a sparse vector of MCHits to accumulate energy in each cube.  But that's a map, you say!  
    //std::map is a more memory-efficient way to implement a sparse vector than just a std::vector with lots of blank 
    //entries.  Think of the RAM needed for ~1e7 MCHits in each event!  
    std::map<GridHits::Triple, GridHits::HitData> hits;
                                                                                                                         
    //Set up to determine whether each TG4HitSegment came from a neutron 
    const auto neutDescendIDs = NeutDescend();
    if(neutDescendIDs.empty()) return false; //If there are no neutron-descneded hits in this event, there is nothing to do.

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    { 
      for(const auto& seg: det.second)
      {
        //Fiducial cut
        const auto start = geo::InLocal(seg.Start.Vect(), mat);
        const auto stop = geo::InLocal(seg.Stop.Vect(), mat);
        double arr[] = {start.X(), start.Y(), start.Z()};
        if(shape->Contains(arr)) 
        {
          fHitAlg.MakeHitData(seg, hits, mat, [&neutDescendIDs](const auto& seg){ return !(neutDescendIDs.count(seg.PrimaryId)); });
        } //If this hit segment is in the fiducial volume
      } //Loop over all hit segments in this sensitive detector 
    } //For each sensitive detector

    //Group hits by whether they passed the neighbor cut.  Then, I can perform another neighbor cut among neutron-caused hits to 
    //weed out hits that are part of neutron-induced tracks that start too close to non-neutron or non-visible hits.  
    std::map<GridHits::Triple, std::list<GridHits::Triple>> passedHits, failedHits;
    for(const auto& pair: hits)
    {
      //const auto out = fHitAlg.MakeHit(pair, mat);
      const auto& hit = pair.second;
      if(hit.Energy > fEMin && hit.Energy > 4.*hit.OtherE) 
      {
        std::list<GridHits::Triple> neutronNeighbors;
        if(Neighbors(pair, hits, fNeighborDist, neutronNeighbors)) 
        {
          //fHits.push_back(out); //Look for adjacent neighbors
          //hits.remove(pair.first); //TODO: Remove this hit from the list of hits to consider when making Neighbors cuts.  
                                     //      What happens with hits I haven't processed yet?  I really need to know the Neighbors() 
                                     //      results for all potential "good neutron hits" before I can make the Neighbors() cut this 
                                     //      way.  
          passedHits[pair.first].insert(passedHits[pair.first].end(), neutronNeighbors.begin(), neutronNeighbors.end());
        }
        else
        {
          failedHits[pair.first].insert(failedHits[pair.first].end(), neutronNeighbors.begin(), neutronNeighbors.end());
        }
      }
    }

    //Remove all hits that are neighbors of a bad hit, including hits that become bad hits in this process
    //While number of passed hits is changing
    size_t prevSize = passedHits.size()+1;
    while(passedHits.size() < prevSize)
    {
      prevSize = passedHits.size();
      for(const auto& hit: failedHits)
      {
        for(const auto& pos: hit.second) 
        {
          auto found = passedHits.find(pos);
          if(found != passedHits.end())
          {
            failedHits[pos] = found->second;
            passedHits.erase(found);
          }
          //::RecursiveRemove(pos, passedHits, failedHits);
        }
        failedHits.erase(hit.first);
      }
    }

    //Save the remaining hits that were caused primarily by ancestors of FS neutrons and were isolated from hits that will not be saved.  
    for(const auto& good: passedHits)
    {
      fHits.push_back(fHitAlg.MakeHit(*(hits.find(good.first)), mat));
    }

    return !(fHits.empty());
  }

  std::set<int> GridNeutronHits::NeutDescend()
  {
    std::set<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto& trajs = fEvent->Trajectories;
    const auto& vertices = fEvent->Primaries;
    for(const auto& vtx: vertices)
    {
      for(const auto& prim: vtx.Particles)
      {
        const auto mom = trajs[prim.TrackId].InitialMomentum;
        if(prim.Name == "neutron" && mom.E()-mom.Mag() > fEMin)
        {
          truth::Descendants(prim.TrackId, trajs, neutDescendIDs);
          neutDescendIDs.insert(prim.TrackId);
        }
        //else std::cout << "Primary named " << prim.Name << " with KE " << mom.E()-mom.Mag() << " is not a FS neutron.\n";
      }
    }
    return neutDescendIDs;
  }

  //TODO: I *could* unwrap these loops at compile-time, but I don't see a good reason to put in that much effort just yet.  
  //TODO: Make Neighbors cut on hits themselves.  Maybe record somewhere where Neighbors() cut failed so that I don't get an awful recursive mess?
  //Loop over HitData in map of all hits and return whether there is a non-neutron hit within nCubes of cand.   
  bool GridNeutronHits::Neighbors(const std::pair<GridHits::Triple, GridHits::HitData>& cand, 
                                  const std::map<GridHits::Triple, GridHits::HitData>& hits, const size_t nCubes, 
                                  std::list<GridHits::Triple>& neutronNeighbors) const
  {
    bool noNeighbors = true;
    auto key = cand.first;
    for(int xOff = -(int)nCubes; xOff < (int)nCubes+1; ++xOff)
    {
      for(int yOff = -(int)nCubes; yOff < (int)nCubes+1; ++yOff)
      {
        for(int zOff = -(int)nCubes; zOff < (int)nCubes+1; ++zOff)
        {
          auto offPos = key;
          offPos.First += xOff;
          offPos.Second += yOff;
          offPos.Third += zOff;
          auto found = hits.find(offPos);
          if(found != hits.end() && found->second.Energy > fEMin)
          {
            if(found->second.Energy > 4.*found->second.OtherE) neutronNeighbors.push_back(offPos);
            else noNeighbors = false;
          }
        }
      }
    }
    return noNeighbors;
  }

  REGISTER_PLUGIN(GridNeutronHits, plgn::Reconstructor);
}
