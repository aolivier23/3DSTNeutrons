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

      //The position of the center of this MCCluster along with its' starting time.
      TLorentzVector Position; 

      //Parmeters to specifc the size of this MCCluster.  For now, MCClusters are XWidth x YWidth x ZWidth boxes.
      double XWidth;
      double YWidth;
      double ZWidth;

      ClassDef(MCCluster, 1);
  };
}

#endif //PERS_MCCLUSTER_H
