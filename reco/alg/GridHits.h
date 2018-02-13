//File: GridHits.h
//Brief: Shared algorithm for turning a collection of TG4HitSegments into a collection of HitData objects.  
//       Each HitData is a proto-MCHit that the user can convert back into an MCHit with a helper function.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//local includes
#include "reco/alg/GeoFunc.h"
#include "persistency/MCHit.h"

//ROOT includes
#include "TGeoBBox.h"

//c++ includes
#include <iostream>

#ifndef RECO_GRIDHITS_H
#define RECO_GRIDHITS_H

namespace pers
{
  class MCHit;
}

namespace reco
{
  class GridHits
  { 
    public: 
      //Elements of the public interface that the user will interact with
      //3 indices combined into one to be used as an index to a std::map.  Indices can be negative so that I can reconstitute positions 
      //more easily.  Proxy for position the way it is used here.
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
        double OtherE;
        double Time;
        std::vector<int> TrackIDs;
        size_t NContrib;
      };

      GridHits(const double width);
      virtual ~GridHits() = default;

      //Public interface
      //Update a map from Triple (= position) to data to make an MCHit.
      template <class FUNC>
      void MakeHitData(const TG4HitSegment& seg, std::map<Triple, HitData>& hitMap, TGeoMatrix* mat, FUNC&& pred) const      
      {
        //Next, add each segment to the hit(s) it enters.  This way, I loop over each segment exactly once.
        //Not actually storing all of the data for an MCHit because Width is the same for all MCHits made by this algorithm 
        //and Position can be reconstituted from a Triple key. 
        TVector3 center(0., 0., 0.);
    
        //Fiducial cut
        const auto start = geo::InLocal(seg.Start.Vect(), mat);
        const auto stop = geo::InLocal(seg.Stop.Vect(), mat);
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
                   
              Triple pos(std::lrint(boxX/fWidth-0.5), std::lrint(boxY/fWidth-0.5), std::lrint(boxZ/fWidth-0.5)); 
              //std::lrint rounds to the nearest integer and casts (hopefully correctly) to that integer.
              auto& hit = hitMap[pos];
              const double dist = LengthInsideBox(seg, boxCenter, mat);
    
              if(dist > 0)
              {
                hit.Time += seg.Start.T(); //Add time for average time at end
                ++hit.NContrib;
                const double length = (seg.Stop.Vect()-seg.Start.Vect()).Mag();
                hit.TrackIDs.push_back(seg.PrimaryId);
                    
                if(dist <= length+1e-5) 
                {
                  hit.Energy += seg.EnergyDeposit*dist/length; //TODO: remove sanity check on distance
                  if(pred(seg)) hit.OtherE += seg.EnergyDeposit*dist/length; //User hook to keep track of energy from "special" segments
                }
                else std::cerr << "Got distance inside hitBox that is greater than this segment's length!  dist is: " << dist << "\nlength is: " << length << "\n";
              }
            } //Loop over x positions on this segment
          } //Loop over y positions on this segment
        } //Loop over z positions on this segment
      }


      //Turn the elements of the map from MakeHitData back into an MCHit.
      pers::MCHit MakeHit(const std::pair<Triple, HitData>& hitData, TGeoMatrix* mat) const;

    protected:
      //Data members
      double fWidth; //The width of the cubes used to make HitData objects and MCHits
      TGeoBBox fHitBox; //The geometry of one MCHit

      //Internal methods
      double LengthInsideBox(const TG4HitSegment& seg, const TVector3& boxCenter, TGeoMatrix* mat) const;
  };
}

#endif //RECO_GRIDHITS_H
