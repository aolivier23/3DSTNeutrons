//File: NoGridNeutronHits.cpp
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
#include "reco/NoGridNeutronHits.h"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"

//c++ includes
#include <set>

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

  double DistFromInside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    begin.GetXYZ(posArr);
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
}

namespace reco
{
  NoGridNeutronHits::NoGridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(100.), fEMin(2.)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("NoGridNeutronHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool NoGridNeutronHits::DoReconstruct() 
  {
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
        if(prim.Name == "neutron" && mom.E()-mom.Mag() > fEMin) Descendants(prim.TrackId, trajs, neutDescendIDs);
      }
    }

    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    {
      std::list<TG4HitSegment> neutSegs;
      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        if(std::find(neutDescendIDs.begin(), neutDescendIDs.end(), seg.PrimaryId) != neutDescendIDs.end()) neutSegs.push_back(seg);
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutSegs.size() > 0)
      {
        const TG4HitSegment seed = *(neutSegs.begin()); //Make sure this is copied and not a reference since I am about to delete it

        //Get geometry information about this detector
        auto mat = ::findMat(fGeo->FindNode(seed.Start.X(), seed.Start.Y(), seed.Start.Z())->GetVolume()->GetName(), *(fGeo->GetTopNode()));
        if(mat == nullptr) throw util::exception("Volume Not Found") << "Could not find transformation matrix for volume " << det.first << "\n";
        
        const auto start = ::InLocal(seed.Start.Vect(), mat);
        const auto stop = ::InLocal(seed.Stop.Vect(), mat);

        pers::MCHit hit;
        hit.Energy = seed.EnergyDeposit;
        hit.TrackIDs.push_back(seed.PrimaryId);
        hit.Width = fWidth;
        hit.Position = (seed.Start+seed.Stop)*0.5;

        TVector3 boxCenter = hit.Position.Vect();

        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.
        //Now, try to never use seed again.  It's probably safe anyway if I copied it correctly, but I don't entirely trust myself. 

        //Look for TG4HitSegments from neutrons that are inside hitBox
        neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat](auto& seg)
                           {
                             //Find out how much of seg's total length is inside this box
                             const double dist = ::DistFromInside(hitBox, ::InLocal(seg.Start.Vect(), mat), 
                                                                  ::InLocal(seg.Stop.Vect(), mat), boxCenter);
                                                                                                                                        
                             if(dist == 0.0) return false; //If this segment is completely outside 
                                                           //the box that contains seed, keep it for later.
                            
                             //Otherwise, add this segments's energy to the MCHit
                             hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                             hit.Energy += seg.EnergyDeposit;
                             //TODO: Energy-weighted position average

                             return true;
                           }); //Looking for segments in the same box
                                                                                                                                        
        if(hit.Energy > fEMin) fHits.push_back(hit);

      } //While neutSegs is non-empty
    } //For each SensDet

    return !(fHits.empty());
  }

  //Put all of the descendants of parent into ids.  This is a recursive function, so be careful how you use it.  
  void NoGridNeutronHits::Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const
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
  REGISTER_PLUGIN(NoGridNeutronHits, plgn::Reconstructor);
}
