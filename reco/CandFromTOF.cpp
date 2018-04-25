//File: CandFromTOF.cpp
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates NeutronCands from the MCClusters it contains.  
//       Tries to stitch together MCClusters that could have come from the same FS neutron based on timing.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/CandFromTOF.h"
#include "reco/alg/GeoFunc.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"

//c++ includes
#include <numeric> //std::accumulate got moved here in c++14

namespace
{
  //Function for how I want to sort Clusters and Seeds by position here
  bool less(const TLorentzVector& first, const TLorentzVector& second, const TLorentzVector& vertPos)
  {
    const auto deltaT = second.T() - first.T();
    if(std::fabs(deltaT) > 0.7) return deltaT < 0;
    return ((first-vertPos).Vect().Mag() < (second-vertPos).Vect().Mag());
  }
}

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<reco::CandFromTOF>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch from which to read pers::MCClusters.  Usually also the name of the cluster-making algorithm.", 
                "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
  }
}

namespace reco
{
  CandFromTOF::CandFromTOF(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fCands(), 
                                                                       fClusters(*(config.Input), (*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fClusterAlgName((*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fTimeRes(config.Options->Get<double>("--time-res")), fPosRes(10.)
  {
    config.Output->Branch("CandFromTOF", &fCands);
  }

  bool CandFromTOF::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?  
    const auto& vertPos = vertex.front().Position; 

    //First, sort clusters by starting time to timing resolution, then distance to vertex  
    //TODO: TTreeReaderArray::Iterator_t is NOT an STL random access iterator, so I can't use std::sort.  Maybe copy to my own STL container for 
    //      sorting?
    //std::sort(fClusters.begin(), fClusters.end(), std::bind(::less, vertPos)); 
                                                  /*[&vertPos](const auto& first, const auto& second)
                                                  {
                                                    const auto deltaT = second.FirstPosition.T() - first.FirstPosition.T();
                                                    if(std::fabs(deltaT) > 0.7) return deltaT < 0;
                                                    return ((first.FirstPosition-vertPos).Vect().Mag() < (second.FirstPosition-vertPos).Vect().Mag());
                                                  });*/

    std::list<pers::NeutronCand> seeds; //Seeds for neutron candidates

    //Next, loop over clusters in order and seed neutron candidates based on clusters that are close enough to each seed to have been caused 
    //by the same FS neutron. 
    for(auto outerClustPos = fClusters.begin(); outerClustPos != fClusters.end(); ++outerClustPos)
    {
      auto& outer = *outerClustPos;

      pers::NeutronCand seed;
      seed.DepositedEnergy = outer.Energy;
      seed.Start = outer.FirstPosition;
      seed.ClusterAlgToIndices[fClusterAlgName].push_back(outerClustPos.fIndex); //TODO: Is each map element sorted?  

      const auto diff = seed.Start-vertPos;
      const auto dist = diff.Vect().Mag();
      const auto deltaT = diff.T();
      const float c = 299.792; //Speed of light = 300 mm/ns
      seed.Beta = dist/deltaT/c;

      const double distUncert = fPosRes/dist; //relative uncertainty in distance for this energy "measurement"
      const double timeUncert = fTimeRes/deltaT; //relative uncertainty in time for this energy "measurement"
      seed.SigmaBeta = seed.Beta*std::sqrt(distUncert*distUncert+timeUncert*timeUncert); //Uncertainty in beta
      //const auto energy = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2

      //Look for other seeds that are close enough to this seed that the same FS neutron could have visited both points.
      seeds.remove_if([this, &seed, &vertPos, &c](const auto& other)
                      {
                        const auto relDiff = seed.Start - other.Start;
                        const auto relBeta = relDiff.Vect().Mag()/relDiff.T()/c;

                        //TODO: Do this for each cluster in other? 
                        const auto& first = ::less(seed.Start, other.Start, vertPos)?seed:other;
                        if(relBeta < first.Beta) //TODO: Take into account energy deposited in other seeds so far. 
                                                 //TODO: Rough parameterization of "invisible" energy loss of neutrons versus distance?
                        {
                          //Merge other into seed
                          //I could write operator +() to do this, but there is no such example in edep-sim.  Keeping algorithms separated 
                          //from persistency objects for now.  
                          seed.DepositedEnergy += other.DepositedEnergy;
                          seed.Start = first.Start;
                          seed.ClusterAlgToIndices.insert(other.ClusterAlgToIndices.begin(), other.ClusterAlgToIndices.end());
                          seed.Beta = first.Beta;
                          seed.SigmaBeta = first.SigmaBeta;

                          return true;
                        }
                        return false;
                      });

      seeds.push_back(seed);
    }

    //Calculate candidate aggregate properties
    for(auto& cand: seeds)
    {
      //Accumulate TrackIDs of Clusters in this candidate
      for(const auto& src: cand.ClusterAlgToIndices)
      {
        for(const auto& index: src.second) 
        {
          const auto& clust = fClusters[index];
          cand.TrackIDs.insert(clust.TrackIDs.begin(), clust.TrackIDs.end());
        }
      }

      //Calculate neutron energy from TOF
      //TODO: Use other clusters to refine energy estimate?
      const double mass = 939.6; //Neutron mass in MeV/c^2
      cand.TOFEnergy = mass/std::sqrt(1.-cand.Beta*cand.Beta); //E = gamma * mc^2

      fCands.push_back(cand);
    }

    return !(fCands.empty());
  }
  REGISTER_PLUGIN(CandFromTOF, plgn::Reconstructor);
}

