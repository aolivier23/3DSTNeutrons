//File: CandFromPDF.cpp
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates NeutronCands from the MCClusters it contains.  
//       Tries to stitch together MCClusters that could have come from the same FS neutron based on timing.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

//EDepNeutrons includes
#include "app/Factory.cpp"
#include "reco/CandFromPDF.h"
#include "reco/alg/GeoFunc.h"

//ROOT includes
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TFile.h"
#include "TH2D.h"

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
  void RegCmdLine<reco::CandFromPDF>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch from which to read pers::MCClusters.  Usually also the name of the cluster-making algorithm.", 
                "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
    opts.AddKey("--PDF-file", "Name of a .root file with a TH2D named BetaVsEDep.  Histogram will be used as a Probability Density Function to "
                              "distnguish clusters from different FS neutrons.", INSTALL_CONFIG_DIR "BetaVsEDep");
  }
}

namespace reco
{
  CandFromPDF::CandFromPDF(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fCands(), 
                                                                       fClusters(*(config.Input), (*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fClusterAlgName((*(config.Options))["--cluster-alg"].c_str()), 
                                                                       fTimeRes(config.Options->Get<double>("--time-res")), fPosRes(10.), 
                                                                       fBetaVsEDep(nullptr)
  {
    config.Output->Branch("CandFromPDF", &fCands);

    //TODO: Throw exception if file or histogram doesn't exist
    auto pdfFile = TFile::Open((*(config.Options))["--PDF-file"].c_str());
    if(!pdfFile) std::cerr << "Failed to find file " << (*(config.Options))["--PDF-file"] << "\n";
    auto hist = (TH2D*)(pdfFile->Get("BetaVsEDep"));
    if(!hist) std::cerr << "Failed to find histogram named BetaVsEDep in file " << (*(config.Options))["--PDF-file"] << "\n";
    fBetaVsEDep.reset((TH2D*)(hist->Clone()));
    if(fBetaVsEDep->Integral() > 1.0) fBetaVsEDep->Scale(1./fBetaVsEDep->Integral()); //Normalize PDF if it's not already normalized
  }

  bool CandFromPDF::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?  
    const auto& vertPos = vertex.front().Position;

    //Produce a time-ordered list of clusters
    std::cout << "fClusters.GetSize() is " << fClusters.GetSize() << "\n";
    std::vector<decltype(fClusters.begin())> clusters;
    for(auto pos = fClusters.begin(); pos != fClusters.end(); ++pos) clusters.push_back(pos); 
    std::sort(clusters.begin(), clusters.end(), [&vertPos](const auto& first, const auto& second) 
                                                { return ::less((*first).FirstPosition, (*second).FirstPosition, vertPos); });

    std::vector<pers::NeutronCand> cands; //NeutronCands formed

    //Physical constants
    const double mass = 939.6;
    const float c = 299.792; //Speed of light = 300 mm/ns

    //Next, loop over clusters in order and look for the best next cluster.  Calculate the likelihood for this combination of clusters 
    //and the end of each loop and stop just before likelihood is smaller than best group of clusters so far.  If this is the biggest group of 
    //clusters that satisfies these criteria, it becomes the best group of clusters so far.  Do this until there are no unused clusters left.
    while(clusters.size() > 0) //TODO: This seems vulnerable to an infinite loop
    {
      std::cout << "Looking for another neutron candidate...\n";
      double bestLikelihood = 0.;
      std::vector<decltype(fClusters.begin())> bestCand;
      for(auto first = clusters.begin(); first < clusters.end(); ++first)
      {
        const auto diff = ((*(*first)).Position - vertPos);
        const double firstProb = fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin((*(*first)).Energy, std::fabs(diff.Vect().Mag()/diff.T()/c)));
        if(firstProb != 0) //Skip if probability to be first cluster is 0
        {
          double currentLikelihood = std::log10(firstProb);
          std::cout << "Initialized currentLikelihood to " << currentLikelihood << " for a cluster with beta = " 
                    << std::fabs(diff.Vect().Mag()/diff.T()/c) 
                    << " and energy deposited = " << (*(*first)).Energy << " from bin " 
                    << fBetaVsEDep->FindBin((*(*first)).Energy, std::fabs(diff.Vect().Mag()/diff.T()/c)) << "\n";

          std::vector<decltype(fClusters.begin())> currentCand = {*first};
          for(auto outer = first; outer < clusters.end();)
          {
            const auto inner = std::max_element(outer+1, clusters.end(), [this, &outer, &c](auto& firstIt, auto& secondIt)
                                                                         {
                                                                           const auto& first = *firstIt; 
                                                                           const auto firstDiff = (first.Position - (*(*outer)).Position);
                                                                           const auto firstBeta = std::fabs(firstDiff.Vect().Mag()/firstDiff.T())/c;

                                                                           const auto& second = *secondIt;
                                                                           const auto secDiff = (second.Position - (*(*outer)).Position);
                                                                           const auto secBeta = std::fabs(secDiff.Vect().Mag()/secDiff.T())/c; 

                                                                           return std::log10(fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin(first.Energy, firstBeta))) < 
                                                                                  std::log10(fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin(second.Energy, secBeta)));
                                                                         });
            if(inner < clusters.end())
            {
              const auto diff = (*(*inner)).Position - (*(*outer)).Position;
              if(fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin((*(*inner)).Energy, diff.Vect().Mag()/diff.T()/c)) == 0) //TODO: double equality?
              {
                std::cout << "Stopping at likelihood of " << currentLikelihood << " because all remaining clusters have probability 0.  "
                          << "currentCand.size() is " << currentCand.size()  << "\n";

                  break; //TODO: One "break" is already too many...
              }

              const double prob = std::log10(fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin((*(*inner)).Energy, diff.Vect().Mag()/diff.T()/c)));
              if(currentLikelihood + prob > bestLikelihood) 
              {
                std::cout << "Stopping at likelihood of " << currentLikelihood << " -> " << currentLikelihood + prob << " because bestLikelihood is " 
                          << bestLikelihood << ".  currentCand.size() is " << currentCand.size()  << "\n";
                break;
              }
              currentLikelihood += prob;
              currentCand.push_back(*inner);
            }
            outer = inner;
          } //Producing a single candidate
          if(bestCand.size() == 0 || currentCand.size() > bestCand.size()) 
          {
            std::cout << "Setting bestLikelihood to " << currentLikelihood << " and bestCand.size() to " << currentCand.size() << "\n";
            bestLikelihood = currentLikelihood;
            bestCand = currentCand;
          }
        } //If probability for this first cluster is not 0
      } //For each possible first cluster

      //Construct a NeutronCand from the best group of clusters found
      std::cout << "bestCand.size() is " << bestCand.size() << ", and bestLikelihood is " << bestLikelihood << "\n";
      pers::NeutronCand cand;
      const auto& first = *(*(bestCand.begin()));
      const auto diff = first.FirstPosition - vertPos;
      const auto dist = diff.Vect().Mag();
      const auto deltaT = diff.T();
      cand.Beta = dist/deltaT/c;
      cand.SigmaBeta = cand.Beta*std::sqrt(10.*10./dist/dist+fTimeRes*fTimeRes/deltaT/deltaT); //10mm position resolution assumed
      cand.Start = first.FirstPosition;
      for(const auto& iter: bestCand)
      {
        cand.DepositedEnergy += (*iter).Energy;
        cand.ClusterAlgToIndices["CandFromPDF"].push_back(iter.fIndex);
      }
      cands.push_back(cand);

      //Remove these clusters from consideration for other candidates
      for(const auto& pos: bestCand) 
      {
        const auto found = std::find(clusters.begin(), clusters.end(), pos);
        if(found != clusters.end()) clusters.erase(found);
      }
    } //While there are clusters remaining

    //Calculate candidate aggregate properties
    for(auto& cand: cands)
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
  REGISTER_PLUGIN(CandFromPDF, plgn::Reconstructor);
}

