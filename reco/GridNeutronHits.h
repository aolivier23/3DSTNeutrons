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
      size_t fNeighborDist; //How far away should I look for interfering neighbors when deciding to keep hits.  

      GridHits fHitAlg; //Algorithm for grouping TG4HitSegments into MCHits 

      //Internal functions
      std::set<int> NeutDescend(); //Should be const, but I think TTreeReaderArray is not const-correct.
      //TODO: The loop in neighbors *could* be unwrapped at compile-time, but I'm not sure it's worth the extreme amount of effort needed.
      bool Neighbors(const std::pair<GridHits::Triple, GridHits::HitData>& cand, const std::map<GridHits::Triple, GridHits::HitData>& hits, const size_t nCubes) const;
  };
}

#endif //RECO_GRIDNEUTRONHITS_H
