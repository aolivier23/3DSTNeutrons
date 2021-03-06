//File: CandRecoStats.cpp
//Brief: CandRecoStats is an Analyzer that plots quantities related to neutron candidates found as MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/CandRecoStats.h"
#include "app/Factory.cpp"
#include "alg/TruthFunc.h"

//util includes
#include "ROOT/Base/TFileSentry.h"

//c++ includes
#include <numeric>

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::CandRecoStats>(opt::CmdLine& opts)
  {
    opts.AddKey("--cand-alg", "Name of the branch of pers::NeutronCands to analyze.  Usually the name of the candidate-making algorithm.", "CandFromTOF");
    opts.AddKey("--E-min", "Minimum energy for a FS neutron to be plotted.  Should match hit-making and cluster-making algorithms.", "1.5");
  }
}*/

namespace ana
{
  CandRecoStats::CandRecoStats(const plgn::Analyzer::Config& config): plgn::Analyzer(config), 
                                                                      fCands(*(config.Reader), 
                                                                             config.Options["--cand-alg"].as<std::string>().c_str()), 
                                                                      fMinEnergy(config.Options["EMin"].as<double>())
  {
    fCandidateEnergy = config.File->make<TH1D>("CandidateEnergy", "Energy Specturm of Neutron Candidates;Energy [MeV];Events",
                                               150, 0, 150);
    fCandPerNeutron = config.File->make<TH1D>("CandPerNeutron", "Number of Candidates per FS Neutron;Neutron Candidates;Neutrons",
                                              20, 0, 20);
    fNCand = config.File->make<TH1I>("NCand", "Number of Candidates per Event;Neutron Candidates;Events", 20, 0, 20);
    fFSNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of FS Neutrons that Produced Candidates;Energy [MeV];Events",
                                               200, 0, 3000);
    fCauseEnergyVsCandEnergy = config.File->make<TH2D>("CauseEnergyVsCandEnergy", 
                                                       "KE of FS Neutrons versus Energies of their Candidates;Candidate Energy [MeV];"
                                                                                  "FS Neutron KE [MeV];FS Neutrons", 100, 0, 100, 200, 0, 200);
    fCandAngleWRTCause = config.File->make<TH1D>("CandAngleWRTCause", "Angle of Candidate w.r.t. InitialMomentum of FS Neutron;"
                                                                      "#Delta#theta_{Cand} [degrees];Candidates", 180, 0., 180.);
    fDistFromVtx = config.File->make<TH1D>("DistFromVertex", "Distance of Closest Candidate to Vertex per FS Neutron;Distance [mm];FS Neutrons", 
                                           350, 0, 5000);

    fCandPerNeutronVsNeutronKE = config.File->make<TH2D>("CandPerNeutronVsNeutronKE", "Number of Candidates for Each FS Neutron "
                                                                                      "versus Neutron KE;Candidates;KE [MeV];FS Neutrons", 
                                                         200, 0, 3000, 10, 0, 10); 
    fAngleVsDistFromVtx = config.File->make<TH2D>("AngleVsDistFromVtx", "Angle Candidate Makes w.r.t. FS Neutron Momentum versus Distance of "
                                                                        "Candidate from vertex;Distance from Vertex [mm];"
                                                                        "Candidate from vertex;#Delta#theta_{Cand} [degrees];Candidates", 
                                                                        350, 0, 5000, 180, 0., 180.);

    fNeutronsPerCand = config.File->make<TH1D>("NeutronsPerCand", "Number of Neutrons per Candidate;Neutrons;Candidates", 
                                               10, 0, 10);

    fNNeutronsResidual = config.File->make<TH1D>("NNeutronsResidual", "Number of True Neutrons - Number of Candidates;N_true - N_cand;Events", 
                                                 10, -5, 5);
   
    fNeutronEResidual = config.File->make<TH1D>("NeutronEResidual", "Relative Energy Lost to Neutrons Not Seen;"
                                                                    "#frac{E_{neutron, true} - E_{neutron, reco}}{E_{neutron, true}};Events", 
                                                150, -3, 3);

    fERecoVsTrue = config.File->make<TH2D>("ERecoVsTrue", "Reconstructed Versus True Total Neutron Energy;E_{true} [MeV];E_{reco} [MeV];Events", 
                                           200, 0, 3000, 200, 0, 3000);

    //fLostNeutronE = config.File->make<TH1D>("LostNeutronE", "Energy of Neutrons that were Grouped into a Candidate;Neutrons;KE [MeV]", 
    //                                        200, 0, 3000);
  }

  void CandRecoStats::DoAnalyze()
  {
    std::map<int, int> TrackIDsToFS; //Map from TrackIDs to FS neutron
    const auto trajs = fEvent->Trajectories;

    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const int pdg = part.GetPDGCode();
        const double KE = part.GetMomentum().E() - part.GetMomentum().Mag();
        const int trackId = part.GetTrackId();
        #else
        const int pdg = part.PDGCode;
        const double KE = part.Momentum.E() - part.Momentum.Mag();
        const int trackId = part.TrackId;
        #endif

        if(pdg == 2112 && KE > fMinEnergy)
        {
          std::set<int> descend;
          truth::Descendants(trackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          descend.insert(trackId);
          for(const auto& id: descend) TrackIDsToFS[id] = trackId; 
        }
      }
    }

    fNCand->Fill(fCands.GetSize()); 

    std::map<int, std::vector<pers::NeutronCand>> FSToCands; //Map of FS neutron to the candidates it is responsible for
    for(const auto& cand: fCands)
    {
      fCandidateEnergy->Fill(cand.DepositedEnergy);
       
      std::set<int> FSIds; //The TrackIDs of FS neutrons responsible for this candidate.  Should almost always be only 1.
      for(const auto& id: cand.TrackIDs) FSIds.insert(TrackIDsToFS[id]); 
      fNeutronsPerCand->Fill(FSIds.size());
      
      auto mostE = std::max_element(FSIds.begin(), FSIds.end(), [&trajs](const auto& first, const auto& second)
                                                                {
                                                                  const auto& firstTraj = trajs[first];
                                                                  const auto& secondTraj = trajs[second];
                                                                  #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                                                  return (firstTraj.GetInitialMomentum().E() - firstTraj.GetInitialMomentum().Mag()) 
                                                                         < (secondTraj.GetInitialMomentum().E() - secondTraj.GetInitialMomentum().Mag());
                                                                  #else
                                                                  return (firstTraj.InitialMomentum.E() - firstTraj.InitialMomentum.Mag())
                                                                         < (secondTraj.InitialMomentum.E() - secondTraj.InitialMomentum.Mag());
                                                                  #endif
                                                                });

      //TODO: Restate this problem.  I think I really want to plot KE for neutrons whose first clusters are not the first cluster in a candidate.
      /*for(const auto& id: FSIds)
      {
        if(*id != mostE) 
        {
          const auto& traj = trajs[*id];
          fLostNeutronE->Fill(traj.InitialMomentum.E() - traj.InitialMomentum.Mag())
        }
      }*/

      if(FSIds.size() > 1) std::cout << "Got " << FSIds.size() << " true neutrons for one candidate in event " << fEvent->EventId << "\n";

      double sumCauseE = 0.; //Sum of energy from all causes of this candidate      
      for(const int neutronID: FSIds) //For each FS neutron TrackID
      {
        FSToCands[neutronID].push_back(cand);
        const auto& neutron = trajs[neutronID];
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const auto& neutronInitMom = neutron.GetInitialMomentum();
        const auto& neutronFirstPos = neutron.Points[0].GetPosition();
        #else
        const auto& neutronInitMom = neutron.InitialMomentum;
        const auto& neutronFirstPos = neutron.Points[0].Position;
        #endif

        sumCauseE += neutronInitMom.E() - neutronInitMom.Mag();

        //Figure out the angle of the candidate w.r.t. the FS neutron's initial momentum
        const auto candVec = (cand.Start-neutronFirstPos).Vect();

        const double angle = std::acos(candVec.Unit().Dot(neutronInitMom.Vect().Unit()))*180./3.1415926535897932384626433832;
        fCandAngleWRTCause->Fill(angle);
        fAngleVsDistFromVtx->Fill(candVec.Mag(), angle);

      }
      fCauseEnergyVsCandEnergy->Fill(cand.DepositedEnergy, sumCauseE); //Candidate energy might sometimes be slightly larger than sum of FS 
    }

    for(const auto& pair: FSToCands) 
    {
      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
      const double KE = trajs[pair.first].GetInitialMomentum().E()-trajs[pair.first].GetInitialMomentum().Mag();
      #else
      const double KE = trajs[pair.first].InitialMomentum.E()-trajs[pair.first].InitialMomentum.Mag();
      #endif 

      fFSNeutronEnergy->Fill(KE);
    }

    const double trueVisNeutronKE = std::accumulate(FSToCands.begin(), FSToCands.end(), 0., [&trajs](auto sum, const auto& pair)
                                                   {
                                                     #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                                     const double initE = trajs[pair.first].GetInitialMomentum().E();
                                                     #else
                                                     const double initE = trajs[pair.first].InitialMomentum.E();
                                                     #endif
                                                     return sum + initE - 939.6;
                                                   });

     const double trueVisNeutronE = std::accumulate(FSToCands.begin(), FSToCands.end(), 0., [&trajs](auto sum, const auto& pair)
                                                   {
                                                     #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                                     const double initE = trajs[pair.first].GetInitialMomentum().E();
                                                     #else
                                                     const double initE = trajs[pair.first].InitialMomentum.E();
                                                     #endif
                                                     return sum + initE;
                                                   });

    const double totalCandKE = std::accumulate(fCands.begin(), fCands.end(), 0., [](auto sum, const auto& cand) { return sum + cand.TOFEnergy - 939.6; });

    const double totalCandE = std::accumulate(fCands.begin(), fCands.end(), 0., [](auto sum, const auto& cand) { return sum + cand.TOFEnergy; });

    fNNeutronsResidual->Fill((int)(FSToCands.size()) - (int)(fCands.GetSize()));
    fNeutronEResidual->Fill((trueVisNeutronE - totalCandE)/trueVisNeutronE);
    fERecoVsTrue->Fill(trueVisNeutronKE, totalCandKE);

    //Now, look for all of the candidates for each FS neutron.  
    for(const auto& FS: FSToCands) 
    {
      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
      const auto& FSPos = trajs[FS.first].Points[0].GetPosition(); 
      const double FSKE = trajs[FS.first].GetInitialMomentum().E()-trajs[FS.first].GetInitialMomentum().Mag();
      #else
      const auto& FSPos = trajs[FS.first].Points[0].Position;
      const double FSKE = trajs[FS.first].InitialMomentum.E()-trajs[FS.first].InitialMomentum.Mag();
      #endif

      const auto& closest = std::min_element(FS.second.begin(), FS.second.end(), [&FSPos](const auto& first, const auto& second)
                                                                                 { 
                                                                                   return   (first.Start-FSPos).Vect().Mag2() 
                                                                                          < (second.Start-FSPos).Vect().Mag2();
                                                                                 });
      fDistFromVtx->Fill((closest->Start-FSPos).Vect().Mag());
      fCandPerNeutron->Fill(FS.second.size());
      if(FS.second.size() > 5) std::cout << "Many-candidate event (" << FS.second.size() << " candidates): " << fEvent->EventId << "\n";
      fCandPerNeutronVsNeutronKE->Fill(FSKE, FS.second.size());
    }
  }

  REGISTER_PLUGIN(CandRecoStats, plgn::Analyzer)
}
