//File: CandFromTOF.cpp
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates NeutronCands from the MCClusters it contains.  
//       Tries to stitch together MCClusters that could have come from the same FS neutron based on timing.
//Author: Andrew Olivier aolivier@ur.rochester.edu

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

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<reco::CandFromTOF>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch from which to read pers::MCClusters.  Usually also the name of the cluster-making algorithm.", 
                "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
  }
}*/

namespace reco
{
  CandFromTOF::CandFromTOF(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fCands(), 
                                                                       fClusters(*(config.Input), 
                                                                                 config.Options["ClusterAlg"].as<std::string>().c_str()), 
                                                                       fClusterAlgName(config.Options["ClusterAlg"].as<std::string>().c_str()), 
                                                                       fTimeRes(config.Options["TimeRes"].as<double>()), fPosRes(10.)
  {
    config.Output->Branch("CandFromTOF", &fCands);
  }

  bool CandFromTOF::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?

	#ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
    const auto& vertPos = vertex.front().GetPosition();
    #else
    const auto& vertPos = vertex.front().Position;
    #endif

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

    //Physical constants
    const double mass = 939.6;
    const float c = 299.792; //Speed of light = 300 mm/ns

    //Next, loop over clusters in order and seed neutron candidates based on clusters that are close enough to each seed to have been caused 
    //by the same FS neutron. 
    for(auto outerClustPos = fClusters.begin(); outerClustPos != fClusters.end(); ++outerClustPos)
    {
      auto& outer = *outerClustPos;
      if(outer.FirstPosition.T() - vertPos.T() > 3.*fTimeRes) //3ns
      //TODO: Tune this cut?
      {

        pers::NeutronCand seed;
        seed.DepositedEnergy = outer.Energy;
        seed.Start = outer.FirstPosition;
        seed.ClusterAlgToIndices[fClusterAlgName].push_back(std::distance(fClusters.begin(), outerClustPos)); 

        const auto diff = seed.Start-vertPos;
        const auto dist = diff.Vect().Mag();
        const auto deltaT = diff.T();
        seed.Beta = dist/deltaT/c;

        const double distUncert = fPosRes/dist; //relative uncertainty in distance for this energy "measurement"
        const double timeUncert = fTimeRes/deltaT; //relative uncertainty in time for this energy "measurement"
        seed.SigmaBeta = seed.Beta*std::sqrt(distUncert*distUncert+timeUncert*timeUncert); //Uncertainty in beta
        //const auto energy = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2

        //Look for other seeds that are close enough to this seed that the same FS neutron could have visited both points.
        seeds.remove_if([this, &seed, &vertPos, &c, &mass](const auto& other)
                        {
                          auto closest = std::min_element(other.ClusterAlgToIndices.find(fClusterAlgName)->second.begin(), 
                                                          other.ClusterAlgToIndices.find(fClusterAlgName)->second.end(), 
                                                          [this, &seed](const auto& first, const auto& second)
                                                          {
                                                            const auto& firstClust = fClusters[first];
                                                            const auto& secondClust = fClusters[second];
                                                            return ::less(firstClust.FirstPosition, secondClust.FirstPosition, seed.Start);
                                                          });

                          const auto relDiff = seed.Start - fClusters[*closest].FirstPosition;
                          const auto relBeta = relDiff.Vect().Mag()/relDiff.T()/c;
        
                          const auto& first = ::less(seed.Start, other.Start, vertPos)?seed:other;
                          //std::cout << "Before correction for energy, beta was " << first.Beta << "\n";

                          //TODO: Clusters need to be time-ordered for this to work 
                          const double predictedE = mass/std::sqrt(1.-first.Beta*first.Beta)
                                                    -std::accumulate(other.ClusterAlgToIndices.find(fClusterAlgName)->second.begin(),
                                                                     other.ClusterAlgToIndices.find(fClusterAlgName)->second.end(), 
                                                                     0., [this, &closest](double sum, const auto& index)
                                                                         {
                                                                           if(fClusters[index].Position.T() > fClusters[*closest].Position.T()) 
                                                                             return sum;
                                                                           return sum+fClusters[index].Energy;
                                                                         });
                          double predictedBeta;
                          //TODO: Running into a problem where predictedBeta is -nan.  Physically, this could be caused by candidates where neutrons 
                          //      "deposit" energy from the nuclei they interact with.  
                          if(predictedE > mass) predictedBeta = std::sqrt(1.-mass*mass/predictedE/predictedE);
                          else predictedBeta = first.Beta;
                          //std::cout << "After correction for energy loss, beta is " << predictedBeta << "\n";
                                                                     
                          if(relBeta - predictedBeta <= first.SigmaBeta && predictedE > mass) //TODO: get uncertainty on predictedBeta instead
                          //TODO: Rough parameterization of "invisible" energy loss of neutrons versus distance?
                          {
                            //Merge other into seed
                            //I could write operator +() to do this, but there is no such example in edep-sim.  Keeping algorithms separated 
                            //from persistency objects for now.  
                            seed.DepositedEnergy += other.DepositedEnergy;
                            seed.Start = first.Start;

                            //merge maps
                            auto& seedMap = seed.ClusterAlgToIndices;
                            const auto& otherMap = other.ClusterAlgToIndices;
                            for(const auto& otherPair: otherMap)
                            {
                              seedMap[otherPair.first].insert(seedMap[otherPair.first].end(), otherPair.second.begin(), otherPair.second.end());
                            }

                            seed.Beta = first.Beta;
                            seed.SigmaBeta = first.SigmaBeta;

                            return true;
                          }
                          return false;
                        });

        seeds.push_back(seed);
      }
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
      cand.TOFEnergy = mass/std::sqrt(1.-cand.Beta*cand.Beta); //E = gamma * mc^2

      fCands.push_back(cand);
    }

    return !(fCands.empty());
  }
  REGISTER_PLUGIN(CandFromTOF, plgn::Reconstructor)
}

