//File: CandFromCluster.cpp
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and create exactly one NeutronCand from each MCCluster it contains.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/CandFromCluster.h"
#include "reco/alg/GeoFunc.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

//c++ includes
#include <numeric> //std::accumulate got moved here in c++14

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<reco::CandFromCluster>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch from which to read pers::MCClusters.  Usually also the name of the cluster-making algorithm.", 
                "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
  }
}

namespace reco
{
  CandFromCluster::CandFromCluster(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fCands(), 
                                                                       fClusters(*(config.Input), (*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fClusterAlgName((*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fTimeRes(config.Options->Get<double>("--time-res")), fPosRes(10.)
  {
    config.Output->Branch("CandFromCluster", &fCands);
  }

  bool CandFromCluster::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?  
    const auto& vertPos = vertex.front().Position; 

    //Physical constants
    const double mass = 939.6;
    const float c = 299.792; //Speed of light = 300 mm/ns

    //Next, loop over clusters in order and seed neutron candidates based on clusters that are close enough to each seed to have been caused 
    //by the same FS neutron. 
    for(auto outerClustPos = fClusters.begin(); outerClustPos != fClusters.end(); ++outerClustPos)
    {
      auto& outer = *outerClustPos;

      pers::NeutronCand seed;
      seed.DepositedEnergy = outer.Energy;
      seed.Start = outer.FirstPosition;
      seed.ClusterAlgToIndices[fClusterAlgName].push_back(outerClustPos.fIndex); 

      const auto diff = seed.Start-vertPos;
      const auto dist = diff.Vect().Mag();
      const auto deltaT = diff.T();
      seed.Beta = dist/deltaT/c;

      const double distUncert = fPosRes/dist; //relative uncertainty in distance for this energy "measurement"
      const double timeUncert = fTimeRes/deltaT; //relative uncertainty in time for this energy "measurement"
      seed.SigmaBeta = seed.Beta*std::sqrt(distUncert*distUncert+timeUncert*timeUncert); //Uncertainty in beta
      //const auto energy = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2

      seed.TrackIDs.insert(outer.TrackIDs.begin(), outer.TrackIDs.end());
      seed.TOFEnergy = mass/std::sqrt(1.-seed.Beta*seed.Beta); //E = gamma * mc^2
      fCands.push_back(seed);
    }

    return !(fCands.empty());
  }
  REGISTER_PLUGIN(CandFromCluster, plgn::Reconstructor);
}

