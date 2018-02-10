//File: NoGridNeutronHits.h
//Brief: A Reconstructor that reads in a TG4Event from edepsim and creates MCHits from energy deposits that were 
//       created by ancestors of FS neutrons. Does not attempt to create MCHits only at grid points.  Instead, each 
//       MCHit is centered around its' seed's midpoint.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/MCHit.h"

#ifndef RECO_NOGRIDNEUTRONHITS_H
#define RECO_NOGRIDNEUTRONHITS_H

namespace reco
{
  class NoGridNeutronHits: public plgn::Reconstructor
  {
    public:
      NoGridNeutronHits(const plgn::Reconstructor::Config& config);
      virtual ~NoGridNeutronHits() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCHits that will be written to tree
      std::vector<pers::MCHit> fHits;
    private:
      //Parameters that I will refer to
      double fWidth; //The width in mm of each dimension of a MCHit. 
      double fEMin; //The energy threshold in MeV for creating an MCHit.  Neutrons 
                    //with less than this amount of KE are not interesting to me.   

      //Internal functions
      void Descendants(const int parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids) const;
  };
}

#endif //RECO_NOGRIDNEUTRONHITS_H
