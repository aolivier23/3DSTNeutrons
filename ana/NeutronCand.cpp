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
                                                                  fMinEnergy(config.Options->Get<double>("--E-min"))
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
  }

  void NeutronCand::DoAnalyze()
  {
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

    //TODO: Surely there is some way to just get the size.  Maybe by subtracting iterators?
    fNCand->Fill(fClusters.GetSize()); 

    std::map<int, std::vector<pers::MCCluster>> FSToCands; //Map of FS neutron to the candidates it is responsible for
    for(const auto& cand: fClusters)
    {
      fCandidateEnergy->Fill(cand.Energy);
       
      std::set<int> FSIds; //The TrackIDs of FS neutrons responsible for this candidate.  Should almost always be only 1.
      for(const auto& id: cand.TrackIDs) FSIds.insert(TrackIDsToFS[id]); 

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
    for(const auto& FS: FSToCands) 
    {
      const auto FSPos = trajs[FS.first].Points[0].Position;
      const auto& closest = std::min_element(FS.second.begin(), FS.second.end(), [&FSPos](const auto& first, const auto& second)
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
