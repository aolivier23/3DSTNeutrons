//File: CCQEChargedFSFilter.h
//Brief: A Reconstructor that reads in a TG4Event from edepsim and returns true if it is a CC0pi event with 0 or 1 protons.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/MCHit.h"
#include "persistency/MCCluster.h"

//ROOT includes
#include "TTreeReaderArray.h"

#ifndef RECO_CCQECHARGEDFSFILTER_H
#define RECO_CCQECHARGEDFSFILTER_H

namespace reco
{
  class CCQEChargedFSFilter: public plgn::Reconstructor
  {
    public:
      CCQEChargedFSFilter(const plgn::Reconstructor::Config& config);
      virtual ~CCQEChargedFSFilter() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.
  };
}

#endif //RECO_CCQECHARGEDFSFILTER_H
