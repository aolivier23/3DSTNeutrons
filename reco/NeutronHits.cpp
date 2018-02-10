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

  std::ostream& operator <<(std::ostream& os, const TVector3& vec)
  {
    os << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
    return os;
  }

  //Return a TVector3 in a different coordinate system
  TVector3 InLocal(const TVector3& pos, TGeoMatrix* mat)
  {
    double master[3] = {}, local[3] = {};
    pos.GetXYZ(master);
    mat->MasterToLocal(master, local);
    return TVector3(local[0], local[1], local[2]);
  }

  //The opposite of InLocal
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
  NeutronHits::NeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(100.), fEMin(2.)
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
      if(traj.Name == "neutron" && mom.E()-mom.Mag() > fEMin) 
      {
        neutDescendIDs.push_back(traj.TrackId);
        Descendants(traj.TrackId, trajs, neutDescendIDs);
      }
    }

    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    //TODO: Fiducial cut
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    { 
      //Get geometry information about this detector
      const std::string fiducial = "volA3DST_PV";
      auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      TVector3 center(); 
      std::list<TG4HitSegment> neutSegs, others;
      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        const auto local = ::InLocal(seg.Start.Vect(), mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr)) //Intentionally not extrapolating to the boundary.  Very reasonable to leave 
                                 //some room before the boundary in a real detector anyway.  
        {
          if(std::find(neutDescendIDs.begin(), neutDescendIDs.end(), seg.PrimaryId) != neutDescendIDs.end()) neutSegs.push_back(seg);
          else others.push_back(seg);
        }
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutSegs.size() > 0)
      {
        TG4HitSegment seed = *(neutSegs.begin());        
        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.

        const auto start = ::InLocal(seed.Start.Vect(), mat);
        const auto stop = ::InLocal(seed.Stop.Vect(), mat);

        //Loop over fWidth-sized cubes that contain some energy from seed.
        //TODO: Count energy from non-neutron-descended particles.  

        //Looping from - to + for ease of update condition, so find out which position is starting point and which is stopping point.  
        //Note that I am potentially looping over some blocks without neutron energy deposits.  I am looping over a bounding box for this line that is 
        //aligned with the detector's axes.  
        const auto xCond = std::minmax({start.X(), stop.X()});
        const auto yCond = std::minmax({start.Y(), stop.Y()});
        const auto zCond = std::minmax({start.Z(), stop.Z()});
        for(double boxX = (std::floor(xCond.first/fWidth)+0.5)*fWidth; boxX < (std::floor(xCond.second/fWidth)+1.0)*fWidth; boxX += fWidth)
        {
          for(double boxY = (std::floor(yCond.first/fWidth)+0.5)*fWidth; boxY < (std::floor(yCond.second/fWidth)+1.0)*fWidth; boxY += fWidth)
          {
            for(double boxZ = (std::floor(zCond.first/fWidth)+0.5)*fWidth; boxZ < (std::floor(zCond.second/fWidth)+1.0)*fWidth; boxZ += fWidth)
            {
              TVector3 boxCenter(boxX, boxY, boxZ);
                
              pers::MCHit hit;
              hit.Energy = 0;
              hit.TrackIDs = std::vector<int>();
              const auto global = ::InGlobal(boxCenter, mat);
              hit.Position = TLorentzVector(global.X(), global.Y(), global.Z(), 0.);
              hit.Width = fWidth;

              size_t nContrib; //The number of segements that contributed to this hit
              neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat, &nContrib](auto& seg)
                                 {
                                   //Find out whether seg is in this box at all.  
                                   if(::DistFromOutside(hitBox, ::InLocal(seg.Start.Vect(), mat), ::InLocal(seg.Stop.Vect(), mat), boxCenter) > 0.0) return false;

                                   //Find out how much of seg's total length is inside this box
                                   const double dist = ::DistFromInside(hitBox, ::InLocal(seg.Start.Vect(), mat), 
                                                                        ::InLocal(seg.Stop.Vect(), mat), boxCenter);
            
                                   ++nContrib;
                                   hit.Position += TLorentzVector(0., 0., 0., seg.Start.T()); //Add time to hit.Position.
                                   const double length = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
                                   hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                                   if(dist > length) //If this segment is entirely inside the same box as seed
                                   {
                                     hit.Energy += seg.EnergyDeposit;
                                     return true;
                                   }

                                   //Set this segment's starting position to where it leaves this box.  
                                   //Prevents double-counting of some of seg's energy.
                                   auto offset = dist*(seg.Stop-seg.Start).Vect().Unit();
                                   seg.Start = seg.Start+TLorentzVector(offset.X(), offset.Y(), offset.Z(), 0.);

                                   hit.Energy += seg.EnergyDeposit*dist/length;
                                   seg.EnergyDeposit = seg.EnergyDeposit*dist/length;
 
                                   return false;
                                 }); //Looking for segments in the same box

              //Make hit.Position.T() the average time
              hit.Position.SetT(hit.Position.T()/nContrib);

              if(hit.Energy > fEMin) 
              {
                //Now, look for energy from segments of non-neutron-descended particles
                double otherE = 0.;
                others.remove_if([&otherE, &hitBox, &boxCenter, mat](auto& seg)
                                 {
                                   if(::DistFromOutside(hitBox, ::InLocal(seg.Start.Vect(), mat), ::InLocal(seg.Stop.Vect(), mat), boxCenter) > 0.0) return false;

                                   //Find out how much of seg's total length is inside this box
                                   const double dist = ::DistFromInside(hitBox, ::InLocal(seg.Start.Vect(), mat), 
                                                                        ::InLocal(seg.Stop.Vect(), mat), boxCenter);

                                   const double length = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
                                   if(dist > length) //If this segment is entirely inside the same box as seed
                                   {
                                     otherE += seg.EnergyDeposit;
                                     return true;
                                   }
 
                                   //Set this segment's starting position to where it leaves this box.
                                   //Prevents double-counting of some of seg's energy.
                                   auto offset = dist*(seg.Stop-seg.Start).Vect().Unit();
                                   seg.Start = seg.Start+TLorentzVector(offset.X(), offset.Y(), offset.Z(), 0.);

                                   otherE += seg.EnergyDeposit*dist/length;
                                   seg.EnergyDeposit = seg.EnergyDeposit*dist/length;

                                   return false;
                                 }); //Looking for segments in the same box

                if(hit.Energy > 3.*otherE) fHits.push_back(hit);
                else std::cout << "Rejected a hit with " << hit.Energy << " MeV because there was " << otherE << " MeV from others.\n";
              } //If hit has more than minimum energy
            } //Loop over box z position
          } //Loop over box y position
        } //Loop over box x position
        //TODO: Remove what is left of seed from the list of neutron hits.  
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
