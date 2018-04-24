//File: MCCluster.h
//Brief: An MCCluster is a group of MCHits.  It accumulates the properties of those MCHits, and it has an energy-weighted center position.
//       In this case, an MCCluster is a box, so it has three widths that are its maximum extent in each direction.     
//       Following the style of TG4 objects from Clark McGrew's edepsim.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include <TObject.h>
#include <TLorentzVector.h>

#ifndef PERS_MCCLUSTER_H
#define PERS_MCCLUSTER_H

namespace pers
{
  class MCCluster: public TObject
  { 
    public:
      MCCluster(void) {}
      virtual ~MCCluster(); //This needs to be defined in the .cpp file so that I can call ClassImp(?)
      
      //Energy deposited in this MCCluster
      double Energy;

      //TrackIDs that contributed to this MCCluster.  
      std::vector<int> TrackIDs; 

      TLorentzVector Position; //The position of the center of this MCCluster along with its' average time.
      TLorentzVector FirstPosition; //The position of the "first" MCHit in this cluster.  Using "first" = "closest to vertex" for now.  

      //Parmeters to specifc the size of this MCCluster.  Farthest points from Position.
      float XWidth;
      float YWidth;
      float ZWidth;      

      //Indices of MCHits in this cluster for each hit algorithm used to make this MCCluster
      //std::map<std::string, std::vector<size_t>> fHitAlgToIndices; //TODO: Turn these indices into hits/links/iterators?  I could do this like what I think ART does and 
                                                                   //      store everything in a Ptr class template that has a unique index, but I'd rather my file format stay 
                                                                   //      as close to edep-sim as possible.  

      ClassDef(MCCluster, 2); //Bumped to version 2 because I have added data members
  };
}

#endif //PERS_MCCLUSTER_H
