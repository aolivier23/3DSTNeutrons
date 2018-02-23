//File: Analyzer.cpp
//An Analyzer reads information from a edepsim TTree, which may have additional objects in it, and 
//       produces other ROOT objects for output.  This base class will provide access to the TG4Event for 
//       derived classes, but derived classes must get their own access to other objects they want to use.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//local includes
#include "ana/Analyzer.h"

//ROOT includes
#include "TTreeReader.h"
#include "TGeoManager.h"

//util includes
#include "ROOT/Base/TFileSentry.h"
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"

namespace plgn
{
  Analyzer::Analyzer(const Config& config): fEvent(*(config.Reader), "Event"), fGeo(nullptr)
  { 
  }

  void Analyzer::Analyze()
  {
    fGeo = gGeoManager; //TODO: Get TGeoManager from the current file instead?  
    DoAnalyze();
  }
}
