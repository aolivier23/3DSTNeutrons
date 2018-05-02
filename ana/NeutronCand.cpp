//File: NeutronCand.cpp
//Brief: NeutronCand is an Analyzer that plots quantities related to neutron candidates found as MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/NeutronCand.h"
#include "app/Factory.cpp"
#include "alg/TruthFunc.h"

//util includes
#include "ROOT/Base/TFileSentry.h"
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

//c++ includes
#include <numeric>

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::NeutronCand>(opt::CmdLine& opts)
  {
    opts.AddKey("--cluster-alg", "Name of the branch of pers::MClusters to analyze.  Usually the name of the cluster-making algorithm.", "MergedClusters");
    opts.AddKey("--E-min", "Minimum energy for a FS neutron to be plotted.  Should match hit-making and cluster-making algorithms.", "1.5");
  }
}

namespace ana
{
  NeutronCand::NeutronCand(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fClusters(*(config.Reader), (*(config.Options))["--cluster-alg"].c_str()), 
                                                                  fMinEnergy(config.Options->Get<double>("--E-min")), fClusterNumber(-314), 
                                                                  fClustersFromEnd(-314), fDeltaAngle(-314), fEDep(-314), fELeft(-314), 
                                                                  fEFromTOF(-314), fDistFromPrev(-314), fDeltaT(-314), fTrueE(-314)
  {
    fCandidateEnergy = config.File->make<TH1D>("CandidateEnergy", "Energy Specturm of Neutron Candidates;Energy [MeV];Events",
                                               150, 0, 150);
    fCandPerNeutron = config.File->make<TH1D>("CandPerNeutron", "Number of Candidates per FS Neutron;Neutron Candidates;Neutrons",
                                              20, 0, 20);
    fNCand = config.File->make<TH1I>("NCand", "Number of Candidates per Event;Neutron Candidates;Events", 20, 0, 20);
    fFSNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of FS Neutrons that Produced Candidates;Energy [MeV];Events",
                                               200, 0, 3000);
    fCauseEnergyVsCandEnergy = config.File->make<TH2D>("CauseEnergyVsCandEnergy", "KE of FS Neutrons versus Energies of their Candidates;Candidate Energy [MeV];"
                                                                                  "FS Neutron KE [MeV];FS Neutrons", 100, 0, 100, 200, 0, 200);
    fCandAngleWRTCause = config.File->make<TH1D>("CandAngleWRTCause", "Angle of Candidate w.r.t. InitialMomentum of FS Neutron;"
                                                                      "#Delta#theta_{Cand} [degrees];Candidates", 180, 0., 180.);
    fDistFromVtx = config.File->make<TH1D>("DistFromVertex", "Distance of Closest Candidate to Vertex per FS Neutron;Distance [mm];FS Neutrons", 350, 0, 5000);

    fCandPerNeutronVsNeutronKE = config.File->make<TH2D>("CandPerNeutronVsNeutronKE", "Number of Candidates for Each FS Neutron "
                                                                                      "versus Neutron KE;Candidates;KE [MeV];FS Neutrons", 
                                                         200, 0, 3000, 10, 0, 10); 
    fAngleVsDistFromVtx = config.File->make<TH2D>("AngleVsDistFromVtx", "Angle Candidate Makes w.r.t. FS Neutron Momentum versus Distance of "
                                                                        "Candidate from vertex;Distance from Vertex [mm];"
                                                                        "Candidate from vertex;#Delta#theta_{Cand} [degrees];Candidates", 
                                                                        350, 0, 5000, 180, 0., 180.);

    fNeutronsPerCand = config.File->make<TH1D>("NeutronsPerCand", "Number of Neutrons per Candidate;Neutrons;Candidates", 
                                               10, 0, 10);
     
    fClusterNumVsEDep = config.File->make<TH2D>("ClusterNumVsEDep", "Position in Time-Ordering of Clusters versus Vluster Energy Deposited "
                                                                    "per True Neutron;Energy Deposited [MeV];Cluster Number;Neutrons", 
                                                150, 0, 150, 20, 0, 20);
  
    //Set up output TTree
    fLikelihoodTree = config.File->make<TTree>("LikelihoodTree", "Variables for Clusters based on their true neutron parent");
    fLikelihoodTree->Branch("ClusterNumber", &fClusterNumber);
    fLikelihoodTree->Branch("ClustersFromEnd", &fClustersFromEnd);
    fLikelihoodTree->Branch("DeltaAngle", &fDeltaAngle);   
    fLikelihoodTree->Branch("EDep", &fEDep);
    fLikelihoodTree->Branch("ELeft", &fELeft);
    fLikelihoodTree->Branch("EFromTOF", &fEFromTOF);
    fLikelihoodTree->Branch("DistFromPrev", &fDistFromPrev);
    fLikelihoodTree->Branch("DeltaT", &fDeltaT);
    fLikelihoodTree->Branch("TrueE", &fTrueE);
  }

  void NeutronCand::DoAnalyze()
  {
    //Reset tree branches to default values
    fClusterNumber = -314;
    fClustersFromEnd = -314;
    fDeltaAngle = -314;
    fEDep = -314;
    fELeft = -314;
    fEFromTOF = -314;
    fDistFromPrev = -314;
    fDeltaT = -314;
    fTrueE = -314;

    std::map<int, int> TrackIDsToFS; //Map from TrackIDs to FS neutron
    const auto trajs = fEvent->Trajectories;

    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        if(part.PDGCode == 2112 && part.Momentum.E() - part.Momentum.Mag() > fMinEnergy)
        {
          std::set<int> descend;
          truth::Descendants(part.TrackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          descend.insert(part.TrackId);
          for(const auto& id: descend) TrackIDsToFS[id] = part.TrackId; 
        }
      }
    }

    fNCand->Fill(fClusters.GetSize()); 

    std::map<int, std::vector<pers::MCCluster>> FSToCands; //Map of FS neutron to the candidates it is responsible for
    for(const auto& cand: fClusters)
    {
      fCandidateEnergy->Fill(cand.Energy);
       
      std::set<int> FSIds; //The TrackIDs of FS neutrons responsible for this candidate.  Should almost always be only 1.
      for(const auto& id: cand.TrackIDs) FSIds.insert(TrackIDsToFS[id]); 
      fNeutronsPerCand->Fill(FSIds.size());

      double sumCauseE = 0.; //Sum of energy from all causes of this candidate      
      for(const int neutronID: FSIds) //For each FS neutron TrackID
      {
        FSToCands[neutronID].push_back(cand);
        const auto& neutron = trajs[neutronID];

        sumCauseE += neutron.InitialMomentum.E() - neutron.InitialMomentum.Mag();
        //fCauseEnergyVsCandEnergy->Fill(cand.Energy, neutron.InitialMomentum.E() - neutron.InitialMomentum.Mag());
        /*if(cand.Energy > neutron.InitialMomentum.E() - neutron.InitialMomentum.Mag())
          std::cout << "Candidate energy, " << cand.Energy << ", is > cause energy: " 
                    << neutron.InitialMomentum.E() - neutron.InitialMomentum.Mag() << ".  This candidate has " 
                    << FSIds.size() << " unique causes.\n";*/ //These are multi-cause events!

        //Figure out the angle of the candidate w.r.t. the FS neutron's initial momentum
        const auto candVec = (cand.Position-neutron.Points[0].Position).Vect();

        const double angle = std::acos(candVec.Unit().Dot(neutron.InitialMomentum.Vect().Unit()))*180./3.1415926535897932384626433832;
        fCandAngleWRTCause->Fill(angle);
        fAngleVsDistFromVtx->Fill(candVec.Mag(), angle);

      }
      fCauseEnergyVsCandEnergy->Fill(cand.Energy, sumCauseE); //Candidate energy might sometimes be slightly larger than sum of FS 
    }

    for(const auto& pair: FSToCands) fFSNeutronEnergy->Fill(trajs[pair.first].InitialMomentum.E()-trajs[pair.first].InitialMomentum.Mag());

    //Now, look for all of the candidates for each FS neutron.  
    for(auto& FS: FSToCands) 
    {
      const auto FSPos = trajs[FS.first].Points[0].Position;

      //Sort candidates first according to time, then according to distance from vertex to within time resolution
      auto& clusters = FS.second;
      std::sort(clusters.begin(), clusters.end(), [&FSPos](auto& first, auto& second) 
                                                  {
                                                    const auto firstDiff = first.Position - FSPos; 
                                                    const auto secondDiff = second.Position - FSPos;
                                                    if(std::fabs((secondDiff - firstDiff).T()) > 0.7) return firstDiff.T() < secondDiff.T(); 
                                                    //TODO: Make 0.7ns timing resolution a parameter
                                                    return firstDiff.Vect().Mag() < secondDiff.Vect().Mag();
                                                  });

      TLorentzVector prevPrevPos = FSPos;
      TLorentzVector prevPos = FSPos;
      if(clusters.size() > 0)
      {
        const auto diff = clusters[0].Position - FSPos;
        const float c = 299.792; //Speed of light = 300 mm/ns
        const auto beta = diff.Vect().Mag()/diff.T()/c;
        const auto mass = 939.56563; //MeV/c^2
        fEFromTOF = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2
      }
      const auto eDepTotal = std::accumulate(clusters.begin(), clusters.end(), 0.f, [](const auto sum, const auto& clust) 
                                                                                    { return sum + clust.Energy;});
      float sumE = 0.;
      for(size_t pos = 0; pos < clusters.size(); ++pos)
      {
        fClusterNumVsEDep->Fill(clusters[pos].Energy, pos);
        sumE += clusters[pos].Energy;
        
        //Fill tree for studying likelihood strategies
        fClusterNumber = pos;
        fClustersFromEnd = clusters.size() - pos;
        fDeltaAngle = (prevPos - prevPrevPos).Vect().Unit().Dot((clusters[pos].Position - prevPos).Vect().Unit());
        fEDep = clusters[pos].Energy;
        fELeft = eDepTotal-sumE;
        fDistFromPrev = (prevPos - clusters[pos].Position).Vect().Mag();
        fDeltaT = (clusters[pos].Position - prevPos).T();
        fTrueE = trajs[FS.first].InitialMomentum.E();
        fLikelihoodTree->Fill();

        prevPrevPos = prevPos;
        prevPos = clusters[pos].Position;
      }

      const auto closest = std::min_element(FS.second.begin(), FS.second.end(), [&FSPos](const auto& first, const auto& second)
                                                                                 { 
                                                                                   return   (first.Position-FSPos).Vect().Mag2() 
                                                                                          < (second.Position-FSPos).Vect().Mag2();
                                                                                 });
      fDistFromVtx->Fill((closest->Position-FSPos).Vect().Mag());
      fCandPerNeutron->Fill(FS.second.size());
      if(FS.second.size() > 5) std::cout << "Many-candidate event (" << FS.second.size() << " candidates): " << fEvent->EventId << "\n";
      fCandPerNeutronVsNeutronKE->Fill(trajs[FS.first].InitialMomentum.E()-trajs[FS.first].InitialMomentum.Mag(), FS.second.size());
    }
  }

  REGISTER_PLUGIN(NeutronCand, plgn::Analyzer);
}
