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
#include "alg/TruthFunc.h"
#include "reco/alg/GeoFunc.h"

//c++ includes
#include <set>

namespace
{
  std::ostream& operator <<(std::ostream& os, const TVector3& vec)
  {
    os << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
    return os;
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
    std::set<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    for(const auto& traj: trajs)
    {
      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
      const auto mom = traj.GetInitialMomentum();
      const auto name = traj.GetName();
      #else
      const auto name = traj.Name;
      const auto mom = traj.InitialMomentum;
      #endif
      if(name == "neutron" && mom.E()-mom.Mag() > fEMin) 
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const int id = traj.GetTrackId();
        #else
        const int id = traj.TrackId;
        #endif
        neutDescendIDs.insert(id);
        truth::Descendants(id, trajs, neutDescendIDs);
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
      auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      TVector3 center(); 
      std::list<TG4HitSegment> neutSegs, others;
      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const auto segStart = seg.GetStart();
        const int segPrim = seg.GetPrimaryId();
        #else
        const auto segStart = seg.Start;
        const int segPrim = seg.PrimaryId;
        #endif

        const auto local = geo::InLocal(segStart.Vect(), mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr)) //Intentionally not extrapolating to the boundary.  Very reasonable to leave 
                                 //some room before the boundary in a real detector anyway.  
        {
          if(neutDescendIDs.count(segPrim)) neutSegs.push_back(seg);
          else others.push_back(seg);
        }
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutSegs.size() > 0)
      {
        TG4HitSegment seed = *(neutSegs.begin());        
        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.

        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const auto seedStart = seed.GetStart();
        const auto seedStop = seed.GetStop();
        #else
        const auto seedStart = seed.Start;
        const auto seedStop = seed.Stop;
        #endif

        const auto start = geo::InLocal(seedStart.Vect(), mat);
        const auto stop = geo::InLocal(seedStop.Vect(), mat);

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
              const auto global = geo::InGlobal(boxCenter, mat);
              hit.Position = TLorentzVector(global.X(), global.Y(), global.Z(), 0.);
              hit.Width = fWidth;

              size_t nContrib; //The number of segements that contributed to this hit
              neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat, &nContrib](auto& seg)
                                 {
                                   #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                   auto& segStart = seg.GetStart();
                                   auto& segStop = seg.GetStop();
                                   const int segPrim = seg.GetPrimaryId();
                                   auto& segEDep = seg.GetEnergyDeposit();
                                   #else
                                   auto& segStart = seg.Start;
                                   auto& segStop = seg.Stop;
                                   const int segPrim = seg.PrimaryId;
                                   auto& segEDep = seg.EnergyDeposit;
                                   #endif

                                   //Find out whether seg is in this box at all.  
                                   if(geo::DistFromOutside(hitBox, geo::InLocal(segStart.Vect(), mat), geo::InLocal(segStop.Vect(), mat), boxCenter) > 0.0) return false;

                                   //Find out how much of seg's total length is inside this box
                                   const double dist = geo::DistFromInside(hitBox, geo::InLocal(segStart.Vect(), mat), 
                                                                        geo::InLocal(segStop.Vect(), mat), boxCenter);
            
                                   ++nContrib;
                                   hit.Position += TLorentzVector(0., 0., 0., segStart.T()); //Add time to hit.Position.
                                   const double length = (segStop.Vect()-segStart.Vect()).Mag();
                                   hit.TrackIDs.push_back(segPrim); //This segment contributed something to this hit
                                   if(dist > length) //If this segment is entirely inside the same box as seed
                                   {
                                     hit.Energy += segEDep;
                                     return true;
                                   }

                                   //Set this segment's starting position to where it leaves this box.  
                                   //Prevents double-counting of some of seg's energy.
                                   auto offset = dist*(segStop-segStart).Vect().Unit();
                                   segStart = segStart+TLorentzVector(offset.X(), offset.Y(), offset.Z(), 0.);

                                   hit.Energy += segEDep*dist/length;
                                   segEDep = segEDep*dist/length;
 
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
                                   #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                   auto& segStart = seg.GetStart();
                                   auto& segStop = seg.GetStop();
                                   auto& segEDep = seg.GetEnergyDeposit();
                                   #else
                                   auto& segStart = seg.Start;
                                   auto& segStop = seg.Stop;
                                   auto& segEDep = seg.EnergyDeposit;
                                   #endif

                                   if(geo::DistFromOutside(hitBox, geo::InLocal(segStart.Vect(), mat), geo::InLocal(segStop.Vect(), mat), boxCenter) > 0.0) return false;

                                   //Find out how much of seg's total length is inside this box
                                   const double dist = geo::DistFromInside(hitBox, geo::InLocal(segStart.Vect(), mat), 
                                                                           geo::InLocal(segStop.Vect(), mat), boxCenter);

                                   const double length = (segStop.Vect()-segStart.Vect()).Mag();
                                   if(dist > length) //If this segment is entirely inside the same box as seed
                                   {
                                     otherE += segEDep;
                                     return true;
                                   }
 
                                   //Set this segment's starting position to where it leaves this box.
                                   //Prevents double-counting of some of seg's energy.
                                   auto offset = dist*(segStop-segStart).Vect().Unit();
                                   segStart = segStart+TLorentzVector(offset.X(), offset.Y(), offset.Z(), 0.);

                                   otherE += segEDep*dist/length;
                                   segEDep = segEDep*dist/length;

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

  REGISTER_PLUGIN(NeutronHits, plgn::Reconstructor);
}
