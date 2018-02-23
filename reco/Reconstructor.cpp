//File: Reconstructor.cpp
//Brief: A Reconstructor is a plugin that looks at the information already in an edepsim TTree and adds to it. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//local includes
#include "reco/Reconstructor.h"

//ROOT includes
#include "TTree.h"
#include "TGeoManager.h"
#include "TTreeReader.h"

//util includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"

namespace plgn
{
  Reconstructor::Reconstructor(const Config& config): fEvent(*(config.Input), "Event"), fGeo(nullptr)
  {
  }

  bool Reconstructor::Reconstruct()
  {
    fGeo = gGeoManager; //TODO: Do I want to retrieve the TGeoManager from the current file instead?  
    return DoReconstruct();
  }
} 
