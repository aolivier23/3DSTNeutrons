//File: Analyzer.h
//Brief: An Analyzer reads information from a edepsim TTree, which may have additional objects in it, and 
//       produces other ROOT objects for output.  This base class will provide access to the TG4Event for 
//       derived classes, but derived classes must get their own access to other objects they want to use.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//edepsim includes
#include "TG4Event.h"

//yaml-cpp includes
#include "yaml-cpp/yaml.h"

//ROOT includes
#include "TTreeReaderValue.h"

class TTreeReader;
class TGeoManager;
class TTree;

namespace util
{
  class TFileSentry;
}

namespace plgn
{
  class Analyzer
  {
    public:
      struct Config //TODO: Provide a constructor so that I can hold reference members?
      {
        util::TFileSentry* File;
        TTreeReader* Reader;
        YAML::Node Options;
      };

      Analyzer(const Config& config);
      virtual ~Analyzer() = default;

      void Analyze(); //Public interface to private implementation

    protected:
      virtual void DoAnalyze() = 0; //Do plotting or other analysis tasks

      TTreeReaderValue<TG4Event> fEvent;
      TGeoManager* fGeo;
  };
}
