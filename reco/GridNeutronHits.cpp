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

  //(Improved?) geometry algorithm:
  //Split space into octants around the interaction vertex.  
  //Determine whether a point is in an octant by asking 3 questions: 
  //Is your (X, Y, Z) < Center.(X, Y, Z)?
  //Assuming boundary group is small, when looking for hit segments in same MCHit, look 
  //only at own octant plus boundary group.  Could make boundary group for each hemisphere if 
  //needed.  
  //
  //If still too many entries, could subdivide again.  Could probably subdivide each octant individually until 
  //few enough hit segments in group.   

  //TODO: Probably terrible object design
  //TODO: If I plan to stick with this class, make std::list<TG4HitSegment> a template parameter and move this to its' own header.

  //Classes for custom geometry sorting algorithm.  
  class ZHemisphere
  {
    public:
      ZHemisphere(const TVector3& center): Center(center) {}
      virtual ~ZHemisphere() = default;
                                                                                                                                                            
      TVector3 Center;
                                                                                                                                                            
      virtual std::pair<TVector3, std::list<TG4HitSegment>>& operator [](const TVector3& pos) = 0;
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
  
      virtual std::pair<TVector3, std::list<TG4HitSegment>>& operator [](const TVector3& pos)
      {
        return (pos.Y() < Center.Y())?(*Minus)[pos]:(*Plus)[pos];
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
  
      std::pair<TVector3, std::list<TG4HitSegment>>& operator [](const TVector3& pos)
      {
        return (pos.X() < Center.X())?Minus[pos]:Plus[pos];
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

      std::pair<TVector3, std::list<TG4HitSegment>> Plus;
      std::pair<TVector3, std::list<TG4HitSegment>> Minus;

      virtual std::pair<TVector3, std::list<TG4HitSegment>>& operator [](const TVector3& pos)
      {
        return (pos.Z() < Center.Z())?Minus:Plus;
      }

      virtual size_t size() const
      {
        return Plus.second.size()+Minus.second.size();
      }

      virtual TG4HitSegment& first()
      {
        if(Plus.second.size() > 0) return *(Plus.second.begin());
        return *(Minus.second.begin());
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
  
      virtual std::pair<TVector3, std::list<TG4HitSegment>>& operator [](const TVector3& pos)
      {
        return (pos.Z() < Center.Z())?Minus[pos]:Plus[pos];
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
  };
}

namespace reco
{
  GridNeutronHits::GridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(100.), fEMin(2.)
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

    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

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
      const size_t subdiv = 7;
      ::XHemisphere neutGeom(center, 1000, subdiv), otherGeom(center, 1000, subdiv); //Split based on the first vertex in this event.  
                                                                                     //1000 is the half-width of the detector.  
                                                                                     //TODO: Get half-width from TGeoBBox.  Maybe even Width TVector3?

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
            neutGeom[(seg.Start+seg.Stop).Vect()*0.5].second.push_back(seg);
          }
          else 
          {
            otherGeom[(seg.Start+seg.Stop).Vect()*0.5].second.push_back(seg);
          }
        }
      }

      //Form MCHits from interesting TG4HitSegments
      while(neutGeom.size() > 0)
      {
        const TG4HitSegment seed = neutGeom.first(); //Make sure this is copied and not a reference since I am about to delete it

        //Remove this and replace the std::lists above where neutGeom and otherGeom are defined to not use crazy geometry algorithm.  
        auto& pair = neutGeom[(seed.Start+seed.Stop).Vect()*0.5];
        const auto center = pair.first; 
        auto& neutSegs = pair.second;
        auto& others = otherGeom[(seed.Start+seed.Stop).Vect()*0.5].second; 

        pers::MCHit hit;
        double timeAvg = 0; //Average time of this MCHit
                            //TODO: Separate hits in time?  

        const double neutE = std::accumulate(neutSegs.begin(), neutSegs.end(), 0., [&hit, &timeAvg](double val, const auto& seg) 
                                                                                   { 
                                                                                     hit.TrackIDs.push_back(seg.PrimaryId);
                                                                                     timeAvg += (seg.Start.T()+seg.Stop.T())*0.5;
                                                                                     return val+seg.EnergyDeposit; 
                                                                                   });
        if(neutE > fEMin)
        {
          const double otherE = std::accumulate(others.begin(), others.end(), 0., [&timeAvg](double val, const auto& seg) 
                                                                                  {
                                                                                    timeAvg += (seg.Start.T()+seg.Stop.T())*0.5;
                                                                                    return val+seg.EnergyDeposit; 
                                                                                  });
  
          if(neutE > 3.*otherE)
          {
            hit.Energy = neutE;
            hit.Width = fWidth;
            hit.Position = TLorentzVector(center.X(), center.y(), center.Z(), timeAvg/(neutSegs.size()+others.size())); 
            fHits.push_back(hit);
          }
        }
        neutSegs.clear();
        others.clear();
      } //While there are still neutron hits to process
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
