//File: GridNeutronHits.h
//Brief: A Reconstructor that reads in a TG4Event from edepsim and creates MCHits from energy deposits that were 
//       created by ancestors of FS neutrons. Does not attempt to create MCHits only at grid points.  Instead, each 
//       MCHit is centered around its' seed's midpoint.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/MCHit.h"
#include "reco/alg/GridHits.h"

#ifndef RECO_GRIDNEUTRONHITS_H
#define RECO_GRIDNEUTRONHITS_H

namespace reco
{
  class GridNeutronHits: public plgn::Reconstructor
  {
    public:
      GridNeutronHits(const plgn::Reconstructor::Config& config);
      virtual ~GridNeutronHits() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCHits that will be written to tree
      std::vector<pers::MCHit> fHits;
    private:
      //Parameters that I will refer to
      double fEMin; //The energy threshold in MeV for creating an MCHit.  Neutrons 
                    //with less than this amount of KE are not interesting to me.  

      GridHits fHitAlg; //Algorithm for grouping TG4HitSegments into MCHits 

      //Internal functions
      std::set<int> NeutDescend() const;
  };
}

#endif //RECO_GRIDNEUTRONHITS_H
