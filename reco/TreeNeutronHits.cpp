//File: TreeNeutronHits.cpp
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
#include "reco/TreeNeutronHits.h"
#include "reco/alg/Octree.cpp"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"
#include "reco/alg/GeoFunc.h"
#include "alg/TruthFunc.h"

//c++ includes
#include <set>
#include <numeric> //TODO: Not needed if not using gcc 6.  std::accumulate should be in <algorithm> instead.

namespace reco
{
  TreeNeutronHits::TreeNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fEMin(2.)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("TreeNeutronHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool TreeNeutronHits::DoReconstruct() 
  {
    std::cout << "Entering TreeNeutronHits::DoReconstruct()\n";
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Remember that I get fEvent and fGeo for free from the base class.  
    //First, figure out which TG4Trajectories are descendants of particles I am interested in.  
    //I am interested in primary neutrons with > 2 MeV KE. 
    std::set<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    const auto vertices = fEvent->Primaries;
    for(const auto& vtx: vertices)
    {
      for(const auto& prim: vtx.Particles)
      {
	#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto primId = prim.GetTrackId();
	const auto mom = trajs[primId].GetInitialMomentum();
        const auto name = prim.GetName();
        #else
        auto primId = prim.TrackId;
	const auto mom = trajs[primId].InitialMomentum;
        const auto name = prim.Name.c_str();
        #endif

        if(strcmp(name, "neutron") == 0 && mom.E()-mom.Mag() > fEMin) 
        {
          truth::Descendants(primId, trajs, neutDescendIDs);
          neutDescendIDs.insert(primId);
        }
        //else std::cout << "Primary named " << prim.Name << " with KE " << mom.E()-mom.Mag() << " is not a FS neutron.\n";
      }
    }
    if(neutDescendIDs.empty()) return false; //If there are no neutron-descneded hits in this event, there is nothing to do.

    std::cout << "Entering loop over SensDets.\n";
    //Next, find all TG4HitSegments that are descended from an interesting FS particle.  
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    {
      //Get geometry information about the detector of interest
      //TODO: This is specfic to the files I am processing right now!
      const std::string fiducial = "volA3DST_PV";
      auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      //Group energy deposits into "subdetectors"
      const auto center = geo::InGlobal(TVector3(0., 0., 0.), mat); //Find the center of this detector
      Octree<pers::MCHit, 6> neutGeom(center, TVector3(1200, 1200, 1000)); //TODO: Get half-width as a TVector3 to form a box instead of a cube 
      std::cout << "Created neutGeom.\n";
      Octree<double, 6> otherGeom(center, TVector3(1200, 1200, 1000));
      std::cout << "Succeeded in creating Octrees.\n";

      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        //Simple fiducial cut.  Should really look at how much of deposit is inside the fiducial volume or something.  
        //Ideally, I'll just get edepsim to do this for me in the future by creating a volume for each scintillator block. 
		#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const auto segStart = seg.GetStart();
		const auto segStop = seg.GetStop();
        const int segPrim = seg.GetPrimaryId();
		auto segEnergy = seg.GetEnergyDeposit();
        #else
        const auto segStart = seg.Start;
		const auto segStop = seg.Stop;
        const int segPrim = seg.PrimaryId;
		auto segEnergy = seg.EnergyDeposit;
        #endif

        const auto local = geo::InLocal((segStart.Vect()+segStop.Vect())*0.5, mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr))
        {
          if(neutDescendIDs.count(segPrim))
          {
            const auto center = (segStart+segStop).Vect()*0.5; //TODO: Really, an Octree that knows to split TG4HitSegments would be ideal here
            auto pair = neutGeom[center];
            if(!pair.second) 
            {
              pair.second = new pers::MCHit();
              pair.second->Position = TLorentzVector(pair.first.X(), pair.first.Y(), pair.first.Z(), segStart.T());
            }

            auto& hit = *(pair.second);
            hit.Energy += segEnergy;
            hit.TrackIDs.push_back(segPrim);
          }
          else 
          {
            auto ptr = otherGeom[(segStart+segStop).Vect()*0.5].second;
            if(!ptr) ptr = new double();
            auto& edep = *ptr;
            edep += segEnergy;
          }
        }
      }

      //Write MCHits that pass cuts to this event
      neutGeom.visitor([&otherGeom, this](const auto ptr)
                       {
                         if(ptr) //Ignore empty (= uninitialized) cells
                         {
                           auto& hit = *ptr; 
                           double otherE = 0.;
                           auto otherPtr = otherGeom[hit.Position.Vect()].second;
                           if(otherPtr) otherE = *otherPtr;
                           if(hit.Energy > 3.*otherE && hit.Energy > fEMin) this->fHits.push_back(hit);
                         }
                       }); 
    } //For each SensDet


    return !(fHits.empty());
  }
  REGISTER_PLUGIN(TreeNeutronHits, plgn::Reconstructor)
}
