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

  /*const TG4Trajectory& Matriarch(const TG4Trajectory& child, const std::vector<TG4Trajectory>& trajs)
  {
    if(child.ParentId == -1) return child;
    return Matriarch(trajs[child.ParentId], trajs);
  } 

  const TG4Trajectory& Matriarch(const TG4HitSegment& seg, const std::vector<TG4Trajectory>& trajs)
  {
    return Matriarch(trajs[seg.PrimaryId], trajs);
  }*/

  //(Improved?) geometry algorithm:
  //Split space into octants around the interaction vertex.  
  //Determine which octant a hit segment is in by "asking" 6 questions of the form:
  //Is your start (X, Y, Z) < center (X, Y, Z)?
  //Is your end (X, Y, Z) < center (X, Y, Z)?
  //Are those the same?  
  //If yes, assign hemisphere.  Do other hemispheres.  
  //If no, assign to general boundary group.  
  //Assuming boundary group is small, when looking for hit segments in same MCHit, look 
  //only at own octant plus boundary group.  Could make boundary group for each hemisphere if 
  //needed.  
  //
  //If still too many entries, could subdivide again.  Could probably subdivide each octant individually until 
  //few enough hit segments in group.   

  //TODO: Probably terrible object design

  //Classes for custom geometry sorting algorithm.  
  class ZHemisphere
  {
    public:
      ZHemisphere(const TVector3& center): Center(center) {}
      virtual ~ZHemisphere() = default;
                                                                                                                                                            
      TVector3 Center;
                                                                                                                                                            
      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit) = 0;
      virtual size_t size() const = 0;
      virtual TG4HitSegment& first() = 0;
  };
   
  class YHemisphere
  {
    public:
      YHemisphere(const TVector3& center, const double width, const int subdiv);
   
      virtual ~YHemisphere() = default;  
                                                                                                                                                            
      TVector3 Center;
      std::unique_ptr<ZHemisphere> Plus;
      std::unique_ptr<ZHemisphere> Minus;
  
      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Y() < Center.Y());
        bool stop = (hit.Stop.Y() < Center.Y());
  
        if(start == stop) return start?(*Minus)[hit]:(*Plus)[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hit.Start.Y() - Center.Y()), stopDiff = std::fabs(hit.Stop.Y() - Center.Y());
        if(startDiff > stopDiff)
        {
          return start?(*Minus)[hit]:(*Plus)[hit];
        }
        return stop?(*Minus)[hit]:(*Plus)[hit];
      }

      size_t size() const
      {
        return Plus->size()+Minus->size();
      }

      TG4HitSegment& first()
      {
        if(Plus->size() > 0) return Plus->first();
        return Minus->first();
      }
  };

  class XHemisphere
  {
    public:
      XHemisphere(const TVector3& center, const double width, const int subdiv): Center(center), Plus(center+TVector3(width/2., 0., 0.), width, subdiv), 
                                                                                 Minus(center-TVector3(width/2., 0., 0.), width, subdiv) {};
      virtual ~XHemisphere() = default;
  
      TVector3 Center;
      YHemisphere Plus;
      YHemisphere Minus;
  
      std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.X() < Center.X());
        bool stop = (hit.Stop.X() < Center.X());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hit.Start.Y() - Center.Y()), stopDiff = std::fabs(hit.Stop.Y() - Center.Y());
        if(startDiff > stopDiff)
        {
          return start?Minus[hit]:Plus[hit];
        }
        return stop?Minus[hit]:Plus[hit];
      }

      size_t size() const
      {
        return Plus.size()+Minus.size();
      }

      TG4HitSegment& first()
      {
        if(Plus.size() > 0) return Plus.first();
        return Minus.first();
      }
  };

  //The final Z hemisphere.  Breaks the recursive operator [] calls.
  class ZEnd: public ZHemisphere
  {
    public:
      ZEnd(const TVector3& center): ZHemisphere(center), Plus(), Minus() {}      
      virtual ~ZEnd() = default;

      std::list<TG4HitSegment> Plus;
      std::list<TG4HitSegment> Minus;

      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Z() < Center.Z());
        bool stop = (hit.Stop.Z() < Center.Z());

        if(start == stop) return start?Minus:Plus; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hit.Start.Y() - Center.Y()), stopDiff = std::fabs(hit.Stop.Y() - Center.Y());
        if(startDiff > stopDiff)
        {
          return start?Minus:Plus;
        }
        return stop?Minus:Plus;
      }

      virtual size_t size() const
      {
        return Plus.size()+Minus.size();
      }

      virtual TG4HitSegment& first()
      {
        if(Plus.size() > 0) return *(Plus.begin());
        return *(Minus.begin());
      }
  };

  //Subdivide again!  Hardcoding two subdivisions max in the constructor for now.
  class  ZMore: public ZHemisphere
  {
    public:
      ZMore(const TVector3& center, const double width, const int subdiv): ZHemisphere(center), Plus(center+TVector3(0., 0., width/2.), width/2., subdiv-1), 
                                                                           Minus(center-TVector3(0., 0., width/2.), width/2., subdiv-1) {} 
      ZMore() = delete;
      virtual ~ZMore() = default;
  
      XHemisphere Plus;
      XHemisphere Minus;
  
      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Z() < Center.Z());
        bool stop = (hit.Stop.Z() < Center.Z());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hit.Start.Y() - Center.Y()), stopDiff = std::fabs(hit.Stop.Y() - Center.Y());
        if(startDiff > stopDiff)
        {
          return start?Minus[hit]:Plus[hit];
        }
        return stop?Minus[hit]:Plus[hit];
      }

      virtual size_t size() const
      {
        return Plus.size()+Minus.size();
      }

      virtual TG4HitSegment& first()
      {
        if(Plus.size() > 0) return Plus.first();
        return Minus.first();
      }
  };

  YHemisphere::YHemisphere(const TVector3& center, const double width, const int subdiv): Center(center), Plus(), Minus()
  {
    if(subdiv > 0)
    {
      Plus.reset(new ZMore(center+TVector3(0., width/2., 0.), width, subdiv));
      Minus.reset(new ZMore(center+TVector3(0., width/2., 0.), width, subdiv));
    }
    else
    {
      Plus.reset(new ZEnd(center+TVector3(0., width/2., 0.)));
      Minus.reset(new ZEnd(center+TVector3(0., width/2., 0.)));
    }
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
    //TODO: ::Matriarch() and Descendants() sometimes give different results!  I am not entirely convinced by the 
    //      results of Descendants(), so trying ::Matriarch() as main method.  
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
    TVector3 center(0., 0., 0.);
    if(vertices.size() > 0) center = vertices[0].Position.Vect();

    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.  
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    {
      //Decide on a threshold of when to use this algorithm?  
      //Only sort energy deposits in this detector. 
      //TODO: Come up with a reasoning for these thresholds
      size_t subdiv = 0;
      if(det.second.size() > 1e2) subdiv = 1;
      if(det.second.size() > 1e3) subdiv = 2;
      if(det.second.size() > 1e4) subdiv = 3;

      //TODO: Only use neutGeom if there are lots of energy deposits from FS neutrons? 
      ::XHemisphere neutGeom(center, 1000, subdiv), otherGeom(center, 1000, subdiv); //Split based on the first vertex in this event.  

      //Get geometry information about the detector of interest
      //TODO: This is specfic to the files I am processing right now!
      const std::string fiducial = "volA3DST_PV";
      auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        //Simple fiducial cut.  Should really look at how much of deposit is inside the fiducial volume or something.  
        //Ideally, I'll just get edepsim to do this for me in the future by creating a volume for each scintillator block.  
        const auto local = ::InLocal((seg.Start.Vect()+seg.Stop.Vect())*0.5, mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr))
        {
          //const auto primary = ::Matriarch(seg, trajs);
          auto found = std::find(neutDescendIDs.begin(), neutDescendIDs.end(), seg.PrimaryId);
          if(found != neutDescendIDs.end())
          {
            neutGeom[seg].push_back(seg);
          }
          else 
          {
            otherGeom[seg].push_back(seg);
          }
        }
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutGeom.size() > 0)
      {
        const TG4HitSegment seed = neutGeom.first(); //Make sure this is copied and not a reference since I am about to delete it

        //Remove this and replace the std::lists above where neutGeom and otherGeom are defined to not use crazy geometry algorithm.  
        auto& neutSegs = neutGeom[seed];
        auto& others = otherGeom[seed]; 

        pers::MCHit hit;
        hit.Energy = seed.EnergyDeposit;
        hit.TrackIDs.push_back(seed.PrimaryId);
        hit.Width = fWidth;
        hit.Position = seed.Start;

        TVector3 boxCenter = hit.Position.Vect();

        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.
        //Now, try to never use seed again.  It's probably safe anyway if I copied it correctly, but I don't entirely trust myself. 

        //Look for TG4HitSegments from neutrons that are inside hitBox
        neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat](auto& seg)
                           {
                             //Find out how much of seg's total length is inside this box
                             const double dist = ::DistFromOutside(hitBox, seg.Start.Vect(),
                                                                   seg.Stop.Vect(), boxCenter);
                                                                                                                                          
                             if(dist > 0.0) return false; //If this segment is completely outside 
                                                          //the box that contains seed, keep it for later.
                              
                             //Otherwise, add this segments's energy to the MCHit
                             hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                             hit.Energy += seg.EnergyDeposit;
                             //std::cout << "Accumulated another neutron seg's energy of " << seg.EnergyDeposit << "\n";
                             //TODO: Energy-weighted position average
  
                             return true;
                           }); //Looking for segments in the same box
    
        if(hit.Energy > fEMin)
        {    
          //Add up the non-neutron-descended energy deposits in this box.
          //At least do this calculation when I know what I'm looking for.  
          double otherE = 0.;
          for(auto otherPtr = others.begin(); otherPtr != others.end() /*&& !(hit.Energy+10. > otherE*3.)*/; ++otherPtr) //Add an extra 10 MeV to be sure this is not 
                                                                                                                         //a floating point precision problem
          {
            const auto& seg = *otherPtr;
            //TODO: Should these be InLocal()?
            //TODO: Use result of DistFromInside to figure out how much energy is contributed.
            if(::DistFromInside(hitBox, seg.Start.Vect(), seg.Stop.Vect(), boxCenter) > 0.0)
            {
              otherE += seg.EnergyDeposit;
            }

            //Simpler algorithm to check for segments that either start or end in this hit.
            /*const auto startDiff = seg.Start.Vect()-boxCenter;
            const auto stopDiff = seg.Stop.Vect()-boxCenter;
            if(std::fabs(startDiff.X()) < fWidth || std::fabs(startDiff.Y()) < fWidth || std::fabs(startDiff.Z()) < fWidth) otherE += seg.EnergyDeposit;
            else if(std::fabs(stopDiff.X()) < fWidth || std::fabs(stopDiff.Y()) < fWidth || std::fabs(stopDiff.Z()) < fWidth) otherE += seg.EnergyDeposit;*/

          }

          if(hit.Energy > otherE*3.) fHits.push_back(hit);
        } //If hit is above energy threshold
      } //While neutSegs is non-empty
    } //For each SensDet
  
    return !(fHits.empty());
  }

  //Put all of the descendants of parent into ids.  This is a recursive function, so be careful how you use it.  
  void NoGridNeutronHits::Descendants(const int parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const
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
