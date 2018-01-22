//File: Reconstructor.h
//Brief: A Reconstructor is a plugin that reads in information from an edepsim TG4Event and adds to the TTree that 
//       event came from.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//edepsim includes
#include "TG4Event.h"

//ROOT includes
#include "TTreeReaderValue.h"

class TTreeReader;
class TTree;
class TGeoManager;

namespace plgn
{
  class Reconstructor
  {
    public:
      struct Config
      {
        TTree* Output;
        TTreeReader* Input;
      };

      Reconstructor(const Config& config);
      virtual ~Reconstructor() = default;

      bool Reconstruct(); //Public interface to private implementation

    protected:
      virtual bool DoReconstruct() = 0; //Look at what is already in the tree and do your own reconstruction.

      TTreeReaderValue<TG4Event> fEvent; //Access to the "current" TG4Event.  You'll just have to trust the driver application.
      TGeoManager* fGeo; //Access to the "current" TGeoManager.  Since I might want to change it at some point, setting it from 
                         //this base class.
  };
}
