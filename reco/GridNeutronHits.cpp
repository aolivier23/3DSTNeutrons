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


namespace plgn
{
  //Set up command line parsing
  template <>
  void RegCmdLine<reco::GridNeutronHits>(opt::CmdLine& opts)
  {
    opts.AddKey("--E-min", "In GridNeutronHits, minimum energy for a hit to be visible.", "1.5");
    opts.AddKey("--cube-size", "In GridNeutronHits, size of cube-shaped subdetectors that will become hits.", "10.");
    opts.AddKey("--neighbor-cut", "Cut that requires no nearby energy deposits.", "2");
    opts.AddKey<opt::Exists>("--after-birks", "Use secondary energy deposit which can be calculated after applying Birks' Law.", "false");
    opts.AddKey("--time-res", "Time resolution of 3DST in ns.  Used to smear times of hits.", "0.7");
  }
}

namespace reco
{
  GridNeutronHits::GridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), 
                                                                               fHitAlg(config.Options->Get<double>("--cube-size"), 
                                                                                       config.Options->Get<bool>("--after-birks"), 
                                                                                       config.Options->Get<double>("--time-res"))
  {
    config.Output->Branch("GridNeutronHits", &fHits);
    
    fEMin = config.Options->Get<double>("--E-min");
    fNeighborDist = config.Options->Get<size_t>("--neighbor-cut");
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

    //Save the hits created if they have a large enough majority of neutron energy
    for(const auto& pair: hits)
    {
      const auto out = fHitAlg.MakeHit(pair, mat);
      const auto& hit = pair.second;
      if(hit.Energy > fEMin && hit.Energy > 4.*hit.OtherE) 
      {
        if(Neighbors(pair, hits, fNeighborDist)) fHits.push_back(out); //Look for adjacent neighbors
      }
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
  bool GridNeutronHits::Neighbors(const std::pair<GridHits::Triple, GridHits::HitData>& cand, 
                                  const std::map<GridHits::Triple, GridHits::HitData>& hits, const size_t nCubes) const
  {
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
          if(found != hits.end() && found->second.Energy > fEMin && !(found->second.Energy > 4.*found->second.OtherE)) return false;
        }
      }
    }
    return true;
  }

  REGISTER_PLUGIN(GridNeutronHits, plgn::Reconstructor);
}
