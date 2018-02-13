//File: GridNeutronHits.cpp
//Brief: A Reconstructor that makes MCHits from energy deposits produced by ancestors of FS neutrons.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Base/exception.h"

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

namespace reco
{
  GridNeutronHits::GridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(10.), fEMin(1.5)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("GridNeutronHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool GridNeutronHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Get geometry information about this detector
    const std::string fiducial = "volA3DST_PV";
    auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();
                                                                                                                         
    //Form MCHits from all remaining hit segments
    //First, create a sparse vector of MCHits to accumulate energy in each cube.  But that's a map, you say!  
    //std::map is a more memory-efficient way to implement a sparse vector than just a std::vector with lots of blank 
    //entries.  Think of the RAM needed for ~1e7 MCHits in each event!  
    std::map<::Triple, ::HitData> hits;
                                                                                                                         
    //Set up to determine whether each TG4HitSegment came from a neutron 
    const auto neutDescendIDs = NeutDescend();
    if(neutDescendIDs.empty()) return false; //If there are no neutron-descneded hits in this event, there is nothing to do.

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    { 
      for(const auto& seg: det.second)
      {
        //Fiducial cut
        const auto start = ::InLocal(seg.Start.Vect(), mat);
        const auto stop = ::InLocal(seg.Stop.Vect(), mat);
        double arr[] = {start.X(), start.Y(), start.Z()};
        const bool isNeut = neutDescendIDs.count(seg.PrimaryId); //Was this energy deposit produced by the descendant of a neutron?
        if(shape->Contains(arr)) 
        {
          fHitAlg.MakeHitData(seg, hits, mat, [&neutDescendIDs](const auto& seg){ return neutDescendIDs.count(seg.PrimarId); });
        } //If this hit segment is in the fiducial volume
      } //Loop over all hit segments in this sensitive detector 
    } //For each sensitive detector

    //Save the hits created if they have a large enough majority of neutron energy
    for(const auto& pair: hits)
    {
      const auto out = fHitAlg.MakeHit(pair, mat)
      if(out.Energy > fEMin && hit.NeutronE > (hit.Energy-hit.NeutronE)*3.) fHits.push_back(out); 
    }

    return !(fHits.empty());
  }

  std::set<int> NeutDescend() const
  {
    std::set<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    const auto vertices = fEvent->Primaries;
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

  REGISTER_PLUGIN(GridNeutronHits, plgn::Reconstructor);
}
