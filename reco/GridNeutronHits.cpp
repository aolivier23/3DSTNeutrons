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
#include "reco/Octree.cpp"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"

//c++ includes
#include <set>
#include <numeric> //TODO: Not needed if not using gcc 6.  std::accumulate should be in <algorithm> instead.

namespace
{
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

  //Return a TVector3 in a different coordinate system
  TVector3 InGlobal(const TVector3& pos, TGeoMatrix* mat)
  {
    double master[3] = {}, local[3] = {};
    pos.GetXYZ(local);
    mat->LocalToMaster(local, master);
    return TVector3(master[0], master[1], master[2]);
  }

  double DistFromInside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    (begin-shapeCenter).GetXYZ(posArr);
    dir.GetXYZ(dirArr);

    //Make sure dir points away from shapeCenter
    if(dir.Dot(shapeCenter-begin) > 0)
    {
      dirArr[0] = -dirArr[0];
      dirArr[1] = -dirArr[1];
      dirArr[2] = -dirArr[2];
    }

    return shape.DistFromInside(posArr, dirArr);
  }

  double DistFromOutside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    (begin-shapeCenter).GetXYZ(posArr);
    dir.GetXYZ(dirArr);

    return shape.DistFromOutside(posArr, dirArr);
  }
}

namespace reco
{
  GridNeutronHits::GridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fEMin(2.)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("GridNeutronHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool GridNeutronHits::DoReconstruct() 
  {
    std::cout << "Entering GridNeutronHits::DoReconstruct()\n";
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Remember that I get fEvent and fGeo for free from the base class.  
    //First, figure out which TG4Trajectories are descendants of particles I am interested in.  
    //I am interested in primary neutrons with > 2 MeV KE. 
    std::vector<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    const auto vertices = fEvent->Primaries;
    for(const auto& vtx: vertices)
    {
      for(const auto& prim: vtx.Particles)
      {
        const auto mom = trajs[prim.TrackId].InitialMomentum;
        if(prim.Name == "neutron" && mom.E()-mom.Mag() > fEMin) 
        {
          Descendants(prim.TrackId, trajs, neutDescendIDs);
          neutDescendIDs.push_back(prim.TrackId);
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
      auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      //Group energy deposits into "subdetectors"
      const auto center = ::InGlobal(TVector3(0., 0., 0.), mat); //Find the center of this detector
      Octree<pers::MCHit, 6> neutGeom(center, TVector3(1200, 1200, 1000)); //TODO: Get half-width as a TVector3 to form a box instead of a cube 
      std::cout << "Created neutGeom.\n";
      Octree<double, 6> otherGeom(center, TVector3(1200, 1200, 1000));
      std::cout << "Succeeded in creating Octrees.\n";

      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        //Simple fiducial cut.  Should really look at how much of deposit is inside the fiducial volume or something.  
        //Ideally, I'll just get edepsim to do this for me in the future by creating a volume for each scintillator block.  
        const auto local = ::InLocal((seg.Start.Vect()+seg.Stop.Vect())*0.5, mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr))
        {
          auto found = std::find(neutDescendIDs.begin(), neutDescendIDs.end(), seg.PrimaryId);
          if(found != neutDescendIDs.end())
          {
            const auto center = (seg.Start+seg.Stop).Vect()*0.5; //TODO: Really, an Octree that knows to split TG4HitSegments would be ideal here
            auto pair = neutGeom[center];
            if(!pair.second) 
            {
              pair.second = new pers::MCHit();
              pair.second->Position = TLorentzVector(pair.first.X(), pair.first.Y(), pair.first.Z(), seg.Start.T());
            }

            auto& hit = *(pair.second);
            hit.Energy += seg.EnergyDeposit;
            hit.TrackIDs.push_back(seg.PrimaryId);
          }
          else 
          {
            auto ptr = otherGeom[(seg.Start+seg.Stop).Vect()*0.5].second;
            if(!ptr) ptr = new double();
            auto& edep = *ptr;
            edep += seg.EnergyDeposit;
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

  //Put all of the descendants of parent into ids.  This is a recursive function, so be careful how you use it.  
  void GridNeutronHits::Descendants(const int parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const
  {
    for(const auto& traj: trajs)
    {
      if(traj.ParentId == parent) 
      {
        const int id = traj.TrackId;
        ids.push_back(id);
        Descendants(id, trajs, ids);
      }
    }
  }
  REGISTER_PLUGIN(GridNeutronHits, plgn::Reconstructor);
}
