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
 
  template <class KEY, class T>
  bool empty(const std::map<KEY, std::vector<T>>& map)
  {
    const auto nonEmpty = std::find_if(map.begin(), map.end(), [](const auto& bin) { return !bin.second.empty(); });
    //if(nonEmpty != map.end()) std::cout << "From ::empty(), map has at least size " << nonEmpty->second.size() << "\n";
    return nonEmpty == map.end();
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
                              "distnguish clusters from different FS neutrons.", INSTALL_CONFIG_DIR "BetaVsEDep.root");
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
    fSmallestProb = fBetaVsEDep->GetEntries(); //1./fBetaVsEDep->GetNbinsX()/fBetaVsEDep->GetNbinsY(); 
  }

  bool CandFromPDF::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?  
    const auto& vertPos = vertex.front().Position;

    //std::cout << "fClusters.GetSize() is " << fClusters.GetSize() << "\n";
    std::vector<pers::NeutronCand> neutrons; //NeutronCands formed

    //Physical constants
    const double mass = 939.6;
    const float c = 299.792; //Speed of light = 300 mm/ns

    //"Combinatorial Kalman Filter" described by https://www.ppd.stfc.ac.uk/Pages/ppd_seminars_170215_talks_dmitry_emeliyanov.pdf
    //First, bin clusters in time resolution-sized bins.  
    std::map<unsigned int, std::vector<decltype(fClusters.begin())>> timeBinnedClusters; //not unordered_map because I want to use these in order
    for(auto cluster = fClusters.begin(); cluster != fClusters.end(); ++cluster)
    {
      timeBinnedClusters[(unsigned int)((*cluster).FirstPosition.T()/fTimeRes)].push_back(cluster);
    }
    const auto nBins = timeBinnedClusters.size();

    for(bool empty = ::empty(timeBinnedClusters); !empty; empty = ::empty(timeBinnedClusters))
    {
      //Next, produce candidates that contain one cluster from each time bin
      //TODO: I am getting an empty cands with at least one element in timeBinnedClusters
      std::vector<std::vector<decltype(fClusters.begin())>> cands;

      const auto first = timeBinnedClusters.begin();
      if(first != timeBinnedClusters.end()) //TODO: This should never happen
      {
        const auto& firstBin = *first;
        //Seed with first bin
        for(const auto& cluster: firstBin.second) cands.push_back({cluster}); 
        cands.emplace(cands.begin()); //Allow for skipping the first layer of clusters

        auto start = timeBinnedClusters.begin();
        ++start;
        for(auto iter = start; iter != timeBinnedClusters.end(); ++iter)
        {
          const auto& bin = *iter;
          for(const auto& cluster: bin.second)
          {
            auto newCands = cands;
            for(auto& cand: cands) 
            {
              auto newCand = cand;
              newCand.push_back(cluster);
              newCands.push_back(newCand); //TODO: Only add cluster to candidate if cluster is within distance light could travel?
            }
            cands = newCands;
          }
        }

        //Remove any completely empty vectors
        cands.erase(cands.begin());
      }
      else std::cerr << "timeBinnedClusters doesn't have any elements!\n";

      //Then, pick the candidate with the best log-likelihood
      const auto best = std::max_element(cands.begin(), cands.end(), [this, &vertPos, &c, &nBins](const auto& first, const auto& second)
                                                                     {
                                                                       double firstLike = 0.;
                                                                       for(const auto& iter: first) 
                                                                       {
                                                                         const auto diff = ((*iter).Position - vertPos);
                                                                         firstLike += std::log10(fBetaVsEDep->GetBinContent(
                                                                         fBetaVsEDep->FindBin((*iter).Energy, 
                                                                                               std::fabs(diff.Vect().Mag()/diff.T()/c))));
                                                                       }
                                                                       //Penalty term for missing bins
                                                                       firstLike += (nBins - first.size())*std::log10(fSmallestProb);

                                                                       double secondLike = 0.;
                                                                       for(const auto& iter: second)
                                                                       {
                                                                         const auto diff = ((*iter).Position - vertPos);
                                                                         secondLike += std::log10(fBetaVsEDep->GetBinContent(
                                                                         fBetaVsEDep->FindBin((*iter).Energy,  
                                                                                               std::fabs(diff.Vect().Mag()/diff.T()/c))));
                                                                       }
                                                                       //Penalty term for missing bins
                                                                       secondLike += (nBins - second.size())*std::log10(fSmallestProb);

                                                                       //return firstLike/first.size() < secondLike/second.size();
                                                                       return firstLike < secondLike; 
                                                                       //TODO: Penalty term instead of likelihood per clusters?  What does this 
                                                                       //      mean?
                                                                     });

      //Construct a NeutronCand from the best group of clusters found
      if(best != cands.end())
      {
        const auto& bestCand = *best;
        std::cout << "bestCand.size() is " << bestCand.size() << "\n";
        pers::NeutronCand neutron;
        const auto& first = *(*(bestCand.begin()));
        const auto diff = first.FirstPosition - vertPos;
        const auto dist = diff.Vect().Mag();
        const auto deltaT = diff.T();
        neutron.Beta = dist/deltaT/c;
        neutron.SigmaBeta = neutron.Beta*std::sqrt(10.*10./dist/dist+fTimeRes*fTimeRes/deltaT/deltaT); //10mm position resolution assumed
        neutron.Start = first.FirstPosition;
        for(const auto& iter: bestCand)
        {
          neutron.DepositedEnergy += (*iter).Energy;
          neutron.ClusterAlgToIndices["CandFromPDF"].push_back(iter.fIndex);
        }
        neutrons.push_back(neutron);

        //Remove these clusters from consideration for other candidates
        for(const auto& pos: bestCand) 
        {
          auto bin = timeBinnedClusters.find((unsigned int)((*pos).FirstPosition.T()/fTimeRes));
          if(bin != timeBinnedClusters.end())
          {
            auto& values = bin->second;
            auto found = std::find(values.begin(), values.end(), pos);
            if(found != values.end()) values.erase(found);
          }
        }
      }
      else std::cerr << "Failed to remove any clusters!  cands.size() is " << cands.size() << ".  timeBinnedClusters.size() is " 
                     << timeBinnedClusters.size() << "\n";
    } //While there are clusters remaining

    //Calculate candidate aggregate properties
    for(auto& neutron: neutrons)
    {
      //Accumulate TrackIDs of Clusters in this candidate
      for(const auto& src: neutron.ClusterAlgToIndices)
      {
        for(const auto& index: src.second) 
        {
          const auto& clust = fClusters[index];
          neutron.TrackIDs.insert(clust.TrackIDs.begin(), clust.TrackIDs.end());
        }
      }

      //Calculate neutron energy from TOF
      //TODO: Use other clusters to refine energy estimate?
      neutron.TOFEnergy = mass/std::sqrt(1.-neutron.Beta*neutron.Beta); //E = gamma * mc^2

      fCands.push_back(neutron);
    }

    return !(fCands.empty());
  }
  REGISTER_PLUGIN(CandFromPDF, plgn::Reconstructor);
}

