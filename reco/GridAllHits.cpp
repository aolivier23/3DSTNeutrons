//File: GridAllHits.cpp
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
#include "reco/GridAllHits.h"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"
#include "reco/alg/GeoFunc.h"

//c++ includes
#include <set>


/*namespace plgn
{
  //Set up command line parsing
  template <>
  void RegCmdLine<reco::GridAllHits>(opt::CmdLine& opts)
  {
    opts.AddKey("--E-min", "In GridNeutronHits, minimum energy for a hit to be visible.", "1.5");
    opts.AddKey("--cube-size", "In GridNeutronHits, size of cube-shaped subdetectors that will become hits.", "10.");
    opts.AddKey("--neighbor-cut", "Cut that requires no nearby energy deposits.", "2");
    opts.AddKey<opt::Exists>("--after-birks", "Use secondary energy deposit which can be calculated after applying Birks' Law.", "false");
    opts.AddKey("--time-res", "Time resolution of 3DST in ns.  Used to smear times of hits.", "0.7");
  }
}*/

namespace reco
{
  GridAllHits::GridAllHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), 
                                                                       fEMin(config.Options["EMin"].as<double>()), 
                                                                       fHitAlg(config.Options["CubeSize"].as<double>(), 
                                                                               config.Options["AfterBirks"].as<bool>(),  
                                                                               config.Options["TimeRes"].as<double>())
  {
    config.Output->Branch("GridAllHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool GridAllHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Get geometry information about this detector
    auto shape = fGeo->FindVolumeFast("volA3DST_PV")->GetShape();
    auto mat = geo::findMat("volA3DST_PV", *(fGeo->GetTopNode()));
    
    //Form MCHits from all remaining hit segments
    //First, create a sparse vector of MCHits to accumulate energy in each cube.  But that's a map, you say!  
    //std::map is a more memory-efficient way to implement a sparse vector than just a std::vector with lots of blank 
    //entries.  Think of the RAM needed for ~1e7 MCHits in each event!  
    std::map<GridHits::Triple, GridHits::HitData> hits; //TODO: Write this interface so I never have to know about this map
                                                                                                                         
    //Next, add each segment to the hit(s) it enters.  This way, I loop over each segment exactly once.
    //Not actually storing all of the data for an MCHit because Width is the same for all MCHits made by this algorithm 
    //and Position can be reconstituted from a Triple key. 

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    { 
      for(const auto& seg: det.second)
      {
        //Fiducial cut
        const auto start = geo::InLocal(seg.Start.Vect(), mat);
        double arr[] = {start.X(), start.Y(), start.Z()};
        if(shape->Contains(arr)) //TODO: Put this back
        {
          fHitAlg.MakeHitData(seg, hits, mat, [](const auto& /*elm*/){ return false; });
        } //If passes fiducial cut
      } //For each segment in this detector
    } //For each detector 

    //Save the hits created
    for(const auto& pair: hits)
    {
      const auto out = fHitAlg.MakeHit(pair, mat);
      if(out.Energy > fEMin) fHits.push_back(out); //TODO: If I were going to make a cut on energy from non-neutrons, this is the place to do it
    }

    return !(fHits.empty());
  }

  REGISTER_PLUGIN(GridAllHits, plgn::Reconstructor);
}
