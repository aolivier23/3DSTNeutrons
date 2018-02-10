//File: GridHits.cpp
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
#include "reco/GridHits.h"
#include "persistency/MCHit.h"
#include "app/Factory.cpp"

//c++ includes
#include <set>

namespace
{
  //3 indices combined into one to be used as an index to a std::map.  Indices can be negative so that I can reconstitute positions 
  //more easily.  
  class Triple
  {
    public:
      Triple(const int first, const int second, const int third): First(first), Second(second), Third(third) {}
      virtual ~Triple() = default;
    
      bool operator <(const Triple& other) const
      {
        if(First != other.First) return First < other.First;
        if(Second != other.Second) return Second < other.Second;
        return Third < other.Third;
      }

      //A 3-vector of indices like (First, Second, Third)
      int First;
      int Second;
      int Third;
  };

  //The data I actually need to save for each MCHit.  The constructor default goes well with std::map::operator[].
  struct HitData
  {
    HitData(): Energy(0.), Time(0.), TrackIDs(), NContrib(0) {}
    virtual ~HitData() = default;
   
    double Energy;
    //double neutronE;
    double Time;
    std::vector<int> TrackIDs;
    size_t NContrib;
  };

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

  //Convert a 3-vector from a global position to the local coordinate system described by mat
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
  GridHits::GridHits(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fHits(), fWidth(10.), fEMin(1e-2)
  {
    //TODO: Rewrite interface to allow configuration?  Maybe pass in opt::CmdLine in constructor, then 
    //      reconfigure from opt::Options after Parse() was called? 

    config.Output->Branch("GridHits", &fHits);
  }

  //Produce MCHits from TG4HitSegments descended from FS neutrons above threshold
  bool GridHits::DoReconstruct() 
  {
    //Get rid of the previous event's MCHits
    fHits.clear();

    //Remember that I get fEvent and fGeo for free from the base class.  
    //Get geometry information for forming MCHits
    TGeoBBox hitBox(fWidth/2., fWidth/2., fWidth/2.);

    //Get geometry information about this detector
    const std::string fiducial = "volA3DST_PV";
    auto mat = findMat(fiducial, *(fGeo->GetTopNode()));
    auto shape = fGeo->FindVolumeFast(fiducial.c_str())->GetShape();
                                                                                                                         
    //Form MCHits from all remaining hit segments
    //First, create a sparse vector of MCHits to accumulate energy in each cube.  But that's a map, you say!  
    //std::map is a more memory-efficient way to implement a sparse vector than just a std::vector with lots of blank 
    //entries.  Think of the RAM needed for ~1e7 MCHits in each event!  
    std::map<::Triple, ::HitData> hits;
                                                                                                                         
    //Next, add each segment to the hit(s) it enters.  This way, I loop over each segment exactly once.
    //Not actually storing all of the data for an MCHit because Width is the same for all MCHits made by this algorithm 
    //and Position can be reconstituted from a Triple key. 
    TVector3 center(0., 0., 0.);

    //Next, find all TG4HitSegments that are descended from an interesting FS particle.   
    //TODO: Fiducial cut
    for(const auto& det: fEvent->SegmentDetectors) //Loop over sensitive detectors
    { 
      for(const auto& seg: det.second)
      {
        //Fiducial cut
        const auto start = ::InLocal(seg.Start.Vect(), mat);
        const auto stop = ::InLocal(seg.Stop.Vect(), mat);
        double arr[] = {start.X(), start.Y(), start.Z()};
        if(shape->Contains(arr)) //TODO: Put this back
        {
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

                //Find out how much of seg's total length is inside this box
                //if(::DistFromOutside(hitBox, ::InLocal(seg.Start.Vect(), mat), ::InLocal(seg.Stop.Vect(), mat), boxCenter) <= 0.0) //TODO: This doesn't work 
                                                                                                                                     //like I expect.  
                //{
                  ::Triple pos(std::lrint(boxX/fWidth-0.5), std::lrint(boxY/fWidth-0.5), std::lrint(boxZ/fWidth-0.5)); 
                  //std::lrint rounds to the nearest integer and casts (hopefully correctly) to that integer.
                  auto& hit = hits[pos];
                  std::cout << "Before adding to it, hit.Energy is " << hit.Energy << "\n";
                  const double dist = ::DistFromInside(hitBox, ::InLocal(seg.Start.Vect(), mat),
                                                       ::InLocal(seg.Stop.Vect(), mat), boxCenter);

                  hit.Time += seg.Start.T(); //Add time for average time at end
                  ++hit.NContrib;
                  const double length = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
                  hit.TrackIDs.push_back(seg.PrimaryId); //This segment contributed something to this hit
                  if(dist < length) hit.Energy += seg.EnergyDeposit*dist/length;
                //} //If seg is inside this box (it doesn't have to be)
              } //Loop over x positions on this segment
            } //Loop over y positions on this segment
          } //Loop over z positions on this segment
        } //If this hit segment is in the fiducial volume
      } //Loop over all hit segments in this sensitive detector 
    } //For each sensitive detector

    //Save the hits created
    for(const auto& pair: hits)
    {
      const auto& key = pair.first;
      const auto& hit = pair.second;
      
      //Reconstitute hit position
      const TVector3 pos(key.First*fWidth+0.5, key.Second*fWidth+0.5, key.Third*fWidth+0.5);
      const auto global = ::InGlobal(pos, mat);

      pers::MCHit out;
      out.Position = TLorentzVector(global.X(), global.Y(), global.Z(), hit.Time/hit.NContrib); //Use average of times of hit segments
      out.Energy = hit.Energy;
      out.Width = fWidth;
      out.TrackIDs = hit.TrackIDs;
      if(out.Energy > fEMin) fHits.push_back(out); //TODO: If I were going to make a cut on energy from non-neutrons, this is the place to do it
    }

    return !(fHits.empty());
  }

  //Put all of the descendants of parent into ids.  This is a recursive function, so be careful how you use it.  
  /*void GridHits::Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const
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
  }*/
  REGISTER_PLUGIN(GridHits, plgn::Reconstructor);
}
