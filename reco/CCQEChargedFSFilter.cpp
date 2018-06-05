//File: CCQEChargedFSFilter.cpp
//Brief: Combines all MCHits that are adjacent to other MCHits into one big cluster.  Then, combines leftover MCHits into clusters that are 5 or fewer hit widths away.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/CCQEChargedFSFilter.h"
#include "reco/alg/GeoFunc.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

//c++ includes

namespace reco
{
  CCQEChargedFSFilter::CCQEChargedFSFilter(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config)
  {
  }

  bool CCQEChargedFSFilter::DoReconstruct()
  {
    //TODO: Write a "Reconstructor" filter for any FS topology.  Probably powered by std::regex.  On a very tight deadline, so 
    //      keeping this as simple as possible for now.  
    std::multiset<std::string> mult;
    for(const auto& vert: fEvent->Primaries)
    {
      for(const auto& part: vert.Particles)
      {
        if(part.Name != "mu-" && part.Name != "mu+" && part.Name != "proton" && part.Name != "neutron") 
        {
          //std::cout << "Found a " << part.Name << ", so returning false from CCQEChargedFSFilter.\n";
          return false;
        }
        mult.insert(part.Name);
      }
    }

    if(mult.count("proton") > 1)
    {
      //std::cout << "Too many protons in CCQEChargedFSFilter, so returning false.\n";
      return false; 
    }
    if(mult.count("mu-") + mult.count("mu+") != 1) 
    {
      //std::cout << mult.count("mu-") << "mu- + " << mult.count("mu+") << " mu+ in CCQEChargedFSFilter, so returning false.\n";
      return false;
    }
    return true;
  }
  REGISTER_PLUGIN(CCQEChargedFSFilter, plgn::Reconstructor);
}

