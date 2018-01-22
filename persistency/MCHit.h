//File: MCHit.h
//Brief: A MCHit is a box with a position that contains all of the energy deposits within its' volume.  
//       Following the style of TG4 objects from Clark McGrew's edepsim.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//edepsim includes
#include "TG4HitSegment.h"

//ROOT includes
#include <TObject.h>
#include <TLorentzVector.h>

#ifndef PERS_MCHIT_H
#define PERS_MCHIT_H

namespace pers
{
  class MCHit: public TObject
  { 
    public:
      MCHit(void) {}
      virtual ~MCHit(); //This needs to be defined in the .cpp file so that I can call ClassImp(?)
      
      //Energy deposited in this MCHit
      double Energy;

      //TrackIDs that contributed to this MCHit.  
      std::vector<int> TrackIDs; 

      //The position of the center of this MCHit along with its' starting time.
      TLorentzVector Position; 

      //Parmeters to specifc the size of this MCHit.  For now, MCHits are Width x Width x Width cubes.
      double Width;

      ClassDef(MCHit, 1);
  };
}

#endif //PERS_MCHIT_H
