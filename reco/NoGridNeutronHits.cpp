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
#include "reco/alg/GeoFunc.h"
#include "alg/TruthFunc.h"

//c++ includes
#include <set>
#include <numeric> //TODO: Not needed if not using gcc 6.  std::accumulate should be in <algorithm> instead.

namespace
{
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
		#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto hitStart = hit.GetStart();
		auto hitStop = hit.GetStop();
        #else
        auto& hitStart = hit.Start;
		auto& hitStop = hit.Stop;
        #endif

        bool start = (hitStart.Y() < Center.Y());
        bool stop = (hitStop.Y() < Center.Y());
  
        if(start == stop) return start?(*Minus)[hit]:(*Plus)[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hitStart.Y() - Center.Y()), stopDiff = std::fabs(hitStop.Y() - Center.Y());
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
		#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto hitStart = hit.GetStart();
		auto hitStop = hit.GetStop();
        #else
        auto& hitStart = hit.Start;
		auto& hitStop = hit.Stop;
        #endif

        bool start = (hitStart.X() < Center.X());
        bool stop = (hitStop.X() < Center.X());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hitStart.Y() - Center.Y()), stopDiff = std::fabs(hitStop.Y() - Center.Y());
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
		#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto hitStart = hit.GetStart();
		auto hitStop = hit.GetStop();
        #else
        auto& hitStart = hit.Start;
		auto& hitStop = hit.Stop;
        #endif

        bool start = (hitStart.Z() < Center.Z());
        bool stop = (hitStop.Z() < Center.Z());

        if(start == stop) return start?Minus:Plus; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hitStart.Y() - Center.Y()), stopDiff = std::fabs(hitStop.Y() - Center.Y());
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
		#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto hitStart = hit.GetStart();
		auto hitStop = hit.GetStop();
        #else
        auto& hitStart = hit.Start;
		auto& hitStop = hit.Stop;
        #endif

        bool start = (hitStart.Z() < Center.Z());
        bool stop = (hitStop.Z() < Center.Z());
  
        if(start == stop) return start?Minus[hit]:Plus[hit]; //A straight line that starts and ends in the same octant cannot visit any other octants

        double startDiff = std::fabs(hitStart.Y() - Center.Y()), stopDiff = std::fabs(hitStop.Y() - Center.Y());
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
  NoGridNeutronHits::NoGridNeutronHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits()
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 
    config.Output->Branch("NoGridNeutronHits", &fHits);

    fEMin = config.Options["EMin"].as<double>();
    fWidth = config.Options["CubeSize"].as<size_t>();
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool NoGridNeutronHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Remember that I get fEvent and fGeo for free from the base class.  
    //First, figure out which TG4Trajectories are descendants of particles I am interested in.  
    //I am interested in primary neutrons with > 2 MeV KE. 
    //TODO: truth::Matriarch() and Descendants() sometimes give different results!  I am not entirely convinced by the 
    //      results of Descendants(), so trying truth::Matriarch() as main method.  
    std::set<int> neutDescendIDs; //TrackIDs of FS neutron descendants
    const auto trajs = fEvent->Trajectories;
    const auto vertices = fEvent->Primaries;
    for(const auto& vtx: vertices)
    {
      for(const auto& prim: vtx.Particles)
      {
	#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const int primId = prim.GetTrackId();
	const auto mom = trajs[primId].GetInitialMomentum();
	const auto name = prim.GetName();
        #else
        const int primId = prim.TrackId;
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
    #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
    auto vertexPos = vertices[0].GetPosition();
    #else
	auto vertexPos = vertices[0].Position;
    #endif	

    TVector3 center(0., 0., 0.);
    if(vertices.size() > 0) center = vertexPos.Vect();

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
      auto mat = geo::findMat(fiducial, *(fGeo->GetTopNode()));
      auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();

      for(const auto& seg: det.second) //Loop over TG4HitSegments in this sensitive detector
      {
        //Simple fiducial cut.  Should really look at how much of deposit is inside the fiducial volume or something.  
        //Ideally, I'll just get edepsim to do this for me in the future by creating a volume for each scintillator block.
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const auto segStart = seg.GetStart();
		const auto segStop = seg.GetStop();
        const int segPrim = seg.GetPrimaryId();
        #else
        const auto segStart = seg.Start;
		const auto segStop = seg.Stop;
        const int segPrim = seg.PrimaryId;
        #endif
  
        const auto local = geo::InLocal((segStart.Vect()+segStop.Vect())*0.5, mat);
        double arr[] = {local.X(), local.Y(), local.Z()};
        if(shape->Contains(arr))
        {
          //const auto primary = truth::Matriarch(seg, trajs);
          if(neutDescendIDs.count(segPrim))
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

        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        auto seedId = seed.GetPrimaryId();
		auto seedEnergy = seed.GetEnergyDeposit();
		const auto seedStart = seed.GetStart();
        #else
        auto seedId = seed.PrimaryId;
		auto seedEnergy = seed.EnergyDeposit;
		const auto seedStart = seed.Start;
        #endif

        pers::MCHit hit;
        hit.Energy = seedEnergy;
        hit.TrackIDs.push_back(seedId);
        hit.Width = fWidth;
        hit.Position = seedStart;

        TVector3 boxCenter = hit.Position.Vect();

        neutSegs.erase(neutSegs.begin()); //remove the seed from the list of neutSegs so that it is not double-counted later.
        //Now, try to never use seed again.  It's probably safe anyway if I copied it correctly, but I don't entirely trust myself. 

        //Look for TG4HitSegments from neutrons that are inside hitBox
        neutSegs.remove_if([&hit, &hitBox, &boxCenter, mat](auto& seg)
                           {
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
                             //Find out how much of seg's total length is inside this box
                             const double dist = geo::DistFromOutside(hitBox, segStart.Vect(),
                                                                   segStop.Vect(), boxCenter);
                                                                                                                                          
                             if(dist > 0.0) return false; //If this segment is completely outside 
                                                          //the box that contains seed, keep it for later.
                              
                             //Otherwise, add this segments's energy to the MCHit
                             hit.TrackIDs.push_back(segPrim); //This segment contributed something to this hit
                             hit.Energy += segEnergy;
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
			#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
		    const auto segStart = seg.GetStart();
			const auto segStop = seg.GetStop();
		    const int segEnergy = seg.GetEnergyDeposit();
		    #else
		    const auto segStart = seg.Start;
			const auto segStop = seg.Stop;
		    const int segEnergy = seg.EnergyDeposit;
		    #endif
            if(geo::DistFromInside(hitBox, segStart.Vect(), segStop.Vect(), boxCenter) > 0.0)
            {
              otherE += segEnergy;
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
  REGISTER_PLUGIN(NoGridNeutronHits, plgn::Reconstructor)
}
