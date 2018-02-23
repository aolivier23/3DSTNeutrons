//File: GridHits.cpp
//Brief: Shared algorithm to turn a collection of TG4HitSegments into MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//edepsim includes
#include "TG4HitSegment.h"

//Include header
#include "reco/alg/GridHits.h"

namespace reco
{
  GridHits::GridHits(const double width): fWidth(width), fHitBox(width/2., width/2., width/2.)
  {
  }
  
  double GridHits::LengthInsideBox(const TG4HitSegment& seg, const TVector3& boxCenter, TGeoMatrix* mat) const
  {
    //TODO: Define start and stop only once for each segment?
    const auto start = geo::InLocal(seg.Start.Vect(), mat);
    const auto stop = geo::InLocal(seg.Stop.Vect(), mat);
    double dist = -1.;
                                                                                                                                       
    //Find out how much of seg's total length is inside this box
    //If this segment starts or ends inside this fHitBox, I need to use DistFromInside() directly.  
    //Otherwise, I need to use DistFromOutside() with both the start and end point to find the distance travelled inside fHitBox.  
    if(geo::Contains(fHitBox, start, boxCenter))
    {
      if(geo::Contains(fHitBox, stop, boxCenter)) //If this segment both starts and stops in fHitBox
      {
        dist = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
      }
      else //Otherwise, this segment leaves fHitBox.  Rely on ROOT to find how far it travels inside.
      {
        dist = geo::DistFromInside(fHitBox, start, stop, boxCenter);
      }
    }
    else if(geo::Contains(fHitBox, stop, boxCenter))
    {
      dist = geo::DistFromInside(fHitBox, stop, start, boxCenter);
    }
    else //This segment passes through fHitBox.  So, use DistFromOutside on both sides of the box and subtract from length to get 
         //distance inside.  If I get a negative distance, then this segment never enters fHitBox to begin with.  
    {
      dist = geo::DistFromOutside(fHitBox, start, stop, boxCenter); //distance from starting point to entering box
      dist += geo::DistFromOutside(fHitBox, stop, start, boxCenter); //distance from stopping point to entering box (in reverse direction)
                                                                                                                                       
      //Computed everything "in place" (almost certainly doesn't matter anyway).  Now, set dist for real. 
      dist = (seg.Stop.Vect()-seg.Start.Vect()).Mag()-dist;
    } //If seg is inside this box (it doesn't have to be)

    return dist;
  }

  pers::MCHit GridHits::MakeHit(const std::pair<Triple, HitData>& hitData, TGeoMatrix* mat) const
  {
    const auto& key = hitData.first;
    const auto& hit = hitData.second;

    //Reconstitute hit position
    const TVector3 pos((key.First+0.5)*fWidth, (key.Second+0.5)*fWidth, (key.Third+0.5)*fWidth);
    const auto global = geo::InGlobal(pos, mat);
      
    pers::MCHit out;
    out.Position = TLorentzVector(global.X(), global.Y(), global.Z(), hit.Time/hit.NContrib); //Use average of times of hit segments
    out.Energy = hit.Energy;
    out.Width = fWidth;
    out.TrackIDs = hit.TrackIDs;
    return out;
  }
}
