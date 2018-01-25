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

  double DistFromOutside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    begin.GetXYZ(posArr);
    dir.GetXYZ(dirArr);

    return shape.DistFromOutside(posArr, dirArr);
  }

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
  //TODO: Way to get Plus, Minus, and Both for Both entries

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
      std::list<TG4HitSegment> Both;
  
      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Y() < Center.Y());
        bool stop = (hit.Stop.Y() < Center.Y());
  
        if(start == stop) return start?(*Minus)[hit]:(*Plus)[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants
        return Both;
      }

      size_t size() const
      {
        return Both.size()+Plus->size()+Minus->size();
      }

      TG4HitSegment& first()
      {
        if(Plus->size() > 0) return Plus->first();
        if(Minus->size() > 0) return Minus->first();
        return *(Both.begin());
      }
  };

  class XHemisphere
  {
    public:
      XHemisphere(const TVector3& center, const double width, const int subdiv): Center(center), Plus(center+TVector3(width/2., 0., 0.), width, subdiv), 
                                                                                 Minus(center-TVector3(width/2., 0., 0.), width, subdiv), 
                                                                                 Both() {};
      virtual ~XHemisphere() = default;
  
      TVector3 Center;
      YHemisphere Plus;
      YHemisphere Minus;
      std::list<TG4HitSegment> Both;
  
      std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.X() < Center.X());
        bool stop = (hit.Stop.X() < Center.X());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants
        return Both;
      }

      size_t size() const
      {
        return Both.size()+Plus.size()+Minus.size();
      }

      TG4HitSegment& first()
      {
        if(Plus.size() > 0) return Plus.first();
        if(Minus.size() > 0) return Minus.first();
        return *(Both.begin());
      }
  };

  //The final Z hemisphere.  Breaks the recursive operator [] calls.
  class ZEnd: public ZHemisphere
  {
    public:
      ZEnd(const TVector3& center): ZHemisphere(center), Plus(), Minus(), Both() {}      
      virtual ~ZEnd() = default;

      std::list<TG4HitSegment> Plus;
      std::list<TG4HitSegment> Minus;
      std::list<TG4HitSegment> Both;

      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Z() < Center.Z());
        bool stop = (hit.Stop.Z() < Center.Z());

        if(start == stop) return start?Minus:Plus; //A straight line that starts and ends in the same octant cannot visit any other octants
        return Both;
      }

      virtual size_t size() const
      {
        return Both.size()+Plus.size()+Minus.size();
      }

      virtual TG4HitSegment& first()
      {
        if(Plus.size() > 0) return *(Plus.begin());
        if(Minus.size() > 0) return *(Minus.begin());
        return *(Both.begin());
      }
  };

  //Subdivide again!  Hardcoding two subdivisions max in the constructor for now.
  class  ZMore: public ZHemisphere
  {
    public:
      ZMore(const TVector3& center, const double width, const int subdiv): ZHemisphere(center), Plus(center+TVector3(0., 0., width/2.), width/2., subdiv-1), 
                                                                           Minus(center-TVector3(0., 0., width/2.), width/2., subdiv-1), Both() {} 
      ZMore() = delete;
      virtual ~ZMore() = default;
  
      XHemisphere Plus;
      XHemisphere Minus;
      std::list<TG4HitSegment> Both;
  
      virtual std::list<TG4HitSegment>& operator [](const TG4HitSegment& hit)
      {
        bool start = (hit.Start.Z() < Center.Z());
        bool stop = (hit.Stop.Z() < Center.Z());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants
        return Both;
      }

      virtual size_t size() const
      {
        return Both.size()+Plus.size()+Minus.size();
      }

      virtual TG4HitSegment& first()
      {
        if(Plus.size() > 0) return Plus.first();
        if(Minus.size() > 0) return Minus.first();
        return *(Both.begin());
      }
  };

  YHemisphere::YHemisphere(const TVector3& center, const double width, const int subdiv): Center(center), Plus(), Minus(), Both()
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
  };
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
    TVector3 center(0., 0., 0.);
    if(vertices.size() > 0) center = vertices[0].Position.Vect();

    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.  
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    {
      //TODO: Decide on a threshold of when to use this algorithm?  
      //      Only sort energy deposits in this detector. 
      ::XHemisphere neutGeom(center, 2400, 1), otherGeom(center, 2400, 1); //Split based on the first vertex in this event.  

      //std::list<TG4HitSegment> neutSegs, others; //Sort HitSegments into "interesting" and "others"
      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        if(std::find(neutDescendIDs.begin(), neutDescendIDs.end(), seg.PrimaryId) != neutDescendIDs.end()) neutGeom[seg].push_back(seg);
        else otherGeom[seg].push_back(seg);
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutGeom.size() > 0)
      {
        const TG4HitSegment seed = neutGeom.first(); //*(neutSegs.begin()); //Make sure this is copied and not a reference since I am about to delete it

        //Remove this and replace the std::lists above where neutGeom and otherGeom are defined to not use crazy geometry algorithm.  
        auto& neutSegs = neutGeom[seed]; //TODO: If seed is in Both, make sure I get everything.  I will probably have to define some kind of 
                                         //      meta-list structure so that I can do the remove operations later in this loop.  In the meantime, 
                                         //      I don't expect Both to be a popular list.  A better solution might be to just make sure a Both 
                                         //      return adds the segment to both lists.  
        auto& others = otherGeom[seed]; //TODO: If seed is in Both, make sure I get everything.  I will probably have to define some kind of 
                                        //      meta-list structure so that I can do the remove operations later in this loop.

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
                             const double dist = ::DistFromOutside(hitBox, ::InLocal(seg.Start.Vect(), mat), 
                                                                   ::InLocal(seg.Stop.Vect(), mat), boxCenter);
                                                                                                                                          
                             if(dist > 0.0) return false; //If this segment is completely outside 
                                                          //the box that contains seed, keep it for later.
                              
                             //Otherwise, add this segments's energy to the MCHit
                             hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                             hit.Energy += seg.EnergyDeposit;
                             //TODO: Energy-weighted position average
  
                             return true;
                           }); //Looking for segments in the same box
        
        //Add up the non-neutron-descended energy deposits in this box.
        //At least do this calculation when I know what I'm looking for.  
        double otherE = 0.;
        for(auto otherPtr = others.begin(); otherPtr != others.end() && 3.*otherE < hit.Energy; ++otherPtr)
        {
          const auto& seg = *otherPtr;
          if(::DistFromOutside(hitBox, ::InLocal(seg.Start.Vect(), mat),
             ::InLocal(seg.Stop.Vect(), mat), boxCenter) <= 0)
          {
            otherE += seg.EnergyDeposit;
          }
        }

        if(hit.Energy > fEMin && (hit.Energy > otherE*3.)) fHits.push_back(hit);
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
