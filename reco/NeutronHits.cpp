//File: NeutronHits.cpp
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
#include "reco/NeutronHits.h"
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
  NeutronHits::NeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(10.), fEMin(2.)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("NeutronHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool NeutronHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Remember that I get fEvent and fGeo for free from the base class.  
    //First, figure out which TG4Trajectories are descendants of particles I am interested in.  
    //I am interested in primary neutrons with > 2 MeV KE.   
    std::vector<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    for(const auto& traj: trajs)
    {
      const auto mom = traj.InitialMomentum;
      if(traj.Name == "neutron" && mom.E()-mom.Mag() > fEMin) Descendants(traj.TrackId, trajs, neutDescendIDs);
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
        TG4HitSegment seed = *(neutSegs.begin()); //Make sure this is copied and not a reference
        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.

        //Get geometry information about this detector
        //TODO: This won't do what I want with subvolumes of the "main" detector around.  Of course, most of this isn't needed 
        //      with subvolumes of the main detector.
        auto mat = ::findMat(fGeo->FindNode(seed.Start.X(), seed.Start.Y(), seed.Start.Z())->GetVolume()->GetName(), *(fGeo->GetTopNode()));
        if(mat == nullptr) throw util::exception("Volume Not Found") << "Could not find transformation matrix for volume " << det.first << "\n";
        
        const auto start = ::InLocal(seed.Start.Vect(), mat);
        const auto stop = ::InLocal(seed.Stop.Vect(), mat);
        const double length = (stop-start).Mag();

        //Loop over fWidth-sized cubes that contain some energy from seed.
        for(double boxX = (((int)(start.X()/fWidth))+0.5)*fWidth; boxX < stop.X(); boxX += fWidth)
        {
          for(double boxY = (((int)(start.Y()/fWidth))+0.5)*fWidth; boxY < stop.Y(); boxY += fWidth)
          {
            for(double boxZ = (((int)(start.Z()/fWidth))+0.5)*fWidth; boxZ < stop.Z(); boxZ += fWidth)
            {
              TVector3 boxCenter(boxX, boxY, boxZ);
                
              pers::MCHit hit;
              hit.Width = fWidth;
              hit.TrackIDs.push_back(seed.PrimaryId);
              hit.Position = TLorentzVector(boxCenter.X(), boxCenter.Y(), boxCenter.Z(), seed.Start.T());

              //Figure out how much of seed's energy was deposited in this box.
              const double dist = ::DistFromInside(hitBox, start, stop, boxCenter);
              hit.Energy = seed.EnergyDeposit*dist/length;

              neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat](auto& seg)
                                 {
                                   //Find out how much of seg's total length is inside this box
                                   const double dist = ::DistFromInside(hitBox, ::InLocal(seg.Start.Vect(), mat), 
                                                                        ::InLocal(seg.Stop.Vect(), mat), boxCenter);
            
                                   if(dist == 0.0) return false; //If this segment is completely outside 
                                                                 //the box that contains seed, keep it for later.

                                   const double length = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
                                   hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                                   if(dist >= length) //If this segment is entirely inside the same box as seed
                                   {
                                     hit.Energy += seg.EnergyDeposit;
                                     return true;
                                   }

                                   //Otherwise, at least part of this segment is not in the same box as seed.  Find out how much is inside.  
                                   hit.Energy += seg.EnergyDeposit*dist/length;
 
                                   //Set this segment's starting position to where it leaves this box.  
                                   //Prevents double-counting of some of seg's energy.
                                   auto offset = dist*(seg.Stop-seg.Start).Vect().Unit();  
                                   seg.Start = seg.Start+TLorentzVector(offset.X(), offset.Y(), offset.Z(), seg.Start.T()); 
                                   //TODO: I think I mangled the time component of seg.Start() here

                                   return false;
                                 }); //Looking for segments in the same box

              if(hit.Energy > fEMin) fHits.push_back(hit);
            } //Loop over box z position
          } //Loop over box y position
        } //Loop over box x position
      } //While there are still netries in neutSegs
    } //For each SensDet

    return !(fHits.empty());
  }

  //Put all of the descendants of parent into ids.  This is a recursive function, so be careful how you use it.  
  void NeutronHits::Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const
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
  REGISTER_PLUGIN(NeutronHits, plgn::Reconstructor);
}
