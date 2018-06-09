//File: CandFromPDF.cpp
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates NeutronCands from the MCClusters it contains.  
//       Tries to stitch together MCClusters that could have come from the same FS neutron based on timing.
//Author: Andrew Olivier aolivier@ur.rochester.edu

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

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<reco::CandFromPDF>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch from which to read pers::MCClusters.  Usually also the name of the cluster-making algorithm.", 
                "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
   
    std::string defaultPath = "";
    const auto value = getenv("EDEPNEUTRONS_CONF_PATH");
    if(value != nullptr) defaultPath = value;
    opts.AddKey("--PDF-file", "Name of a .root file with a TH2D named BetaVsEDep.  Histogram will be used as a Probability Density Function to "
                              "distnguish clusters from different FS neutrons.", defaultPath+"BetaVsEDep.root");
  }
}*/

namespace reco
{
  CandFromPDF::CandFromPDF(const plgn::Reconstructor::Config& config): plgn::Reconstructor(config), fCands(), 
                                                                       fClusters(*(config.Input), 
                                                                                 config.Options["ClusterAlg"].as<std::string>().c_str()), 
                                                                       fClusterAlgName(config.Options["ClusterAlg"].as<std::string>().c_str()), 
                                                                       fTimeRes(config.Options["TimeRes"].as<double>()), fPosRes(10.), 
                                                                       fBetaVsEDep(nullptr)
  {
    config.Output->Branch("CandFromPDF", &fCands);

    //TODO: Throw exception if file or histogram doesn't exist
    const auto fileName = config.Options["PDFFile"].as<std::string>(); //TODO: Check default path
    auto pdfFile = TFile::Open(fileName.c_str());
    if(!pdfFile) std::cerr << "Failed to find file " << fileName << "\n";
    auto hist = (TH2D*)(pdfFile->Get("BetaVsEDep"));
    if(!hist) std::cerr << "Failed to find histogram named BetaVsEDep in file " << fileName << "\n";
    fBetaVsEDep.reset((TH2D*)(hist->Clone()));
    if(fBetaVsEDep->Integral() > 1.0) fBetaVsEDep->Scale(1./fBetaVsEDep->Integral()); //Normalize PDF if it's not already normalized
    fPenaltyTerm = 8.5/fBetaVsEDep->GetEntries(); //2.e-20/fBetaVsEDep->GetNbinsX()/fBetaVsEDep->GetNbinsY(); //1./fBetaVsEDep->GetEntries();  
  }

  bool CandFromPDF::DoReconstruct()
  {
    fCands.clear(); //Clear out the old clusters from last time!

    const auto& vertex = fEvent->Primaries; //TODO: What to do when there are multiple vertices?  
   
    #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
    const auto& vertPos = vertex.front().GetPosition();
    #else
    const auto& vertPos = vertex.front().Position;
    #endif

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
     
    std::cout << "fClusters.GetSize() is " << fClusters.GetSize() << "\n";

    for(bool empty = ::empty(timeBinnedClusters); !empty; empty = ::empty(timeBinnedClusters))
    {
      //Next, produce candidates that contain one cluster from each time bin.  To keep memory usage reasonable, only hold one 
      //candidate in memory at a time.  
      
      //Create a vector of positions in bins.  Represents how the current candidate relates to other possible combinations of 
      //clusters.  Going to the next combination of clusters works somewhat like how I might build a binary adding machine.  
      //Skip bins whose positions are set to the end() iterator.  Also create a vector of end iterators to make loop ending condition simpler.
      std::vector<decltype(timeBinnedClusters.begin()->second.begin())> currentCand, begin, end;
      for(auto& bin: timeBinnedClusters) 
      {
        currentCand.push_back(bin.second.begin()); //Set vector to the first possible combination of clusters
        begin.push_back(bin.second.begin());
        end.push_back(bin.second.end());
      }
       
      //Since I'm looking for the largest likelihood value, keep track of the best group of clusters as well as the best likelihood so far. 
      auto best = currentCand; //TODO: Is this a reasonable default?
      double bestLikelihood = -std::numeric_limits<double>::max(); //TODO: set this to the likelihood of the first candidate instead?

      while(currentCand != end)
      {
        double likelihood = 0.;
        for(size_t pos = 0; pos < currentCand.size() && likelihood > bestLikelihood; ++pos)
        {
          const auto iter = currentCand[pos];
          if(iter == end[pos]) likelihood += std::log10(fPenaltyTerm); //Penalty term for missing bins
          else
          {
            //Otherwise, add the log-likelihood for this bin to the total for this candidate
            const auto diff = (*(*iter)).FirstPosition - vertPos;
            likelihood += std::log10(fBetaVsEDep->GetBinContent(fBetaVsEDep->FindBin(1./std::sqrt((*(*iter)).Energy),
                                     exp(std::fabs(diff.Vect().Mag()/diff.T()/c)))));
          }
        }

        //If this is a better match than the current best match, update the best match
        if(likelihood > bestLikelihood)
        {
          std::cout << "Selecting new best candidate with size " << currentCand.size() 
                    << " and likelihood " << bestLikelihood << " -> " << likelihood << "\n";
          best = currentCand;
          bestLikelihood = likelihood;
        }

        //Update currentCand.  Works somewhat like using several bits to build an adding machine.
        for(int pos = currentCand.size()-1; pos > -1;)
        {
          if(currentCand[pos] == end[pos]) //If this iterator was already end() BEFORE I updated anything, I've already skipped it's bin, 
                                           //so set it back to begin (as if I were handling arbitrary-base multi-digit addition) and update the 
                                           //next iterator by continuing the loop. 
          {
            currentCand[pos] = begin[pos];
            --pos;
          }
          else //If I am not done with looping over this bin, update it.  It may become the un-dereferenceable end() iterator, and my loop should 
               //be ready to handle that by skipping this bin.  
          {
            ++currentCand[pos];
            break;
          }
        }
      }

      //Construct bestCand by skipping end() placeholders in best
      std::vector<decltype(fClusters.begin())> bestCand;
      for(size_t pos = 0; pos < best.size(); ++pos)
      {
        if(best[pos] != end[pos]) bestCand.push_back(*(best[pos]));
      }

      //Construct a NeutronCand from the best group of clusters found
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
        neutron.ClusterAlgToIndices["CandFromPDF"].push_back(std::distance(fClusters.begin(), iter));
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
  REGISTER_PLUGIN(CandFromPDF, plgn::Reconstructor)
}

