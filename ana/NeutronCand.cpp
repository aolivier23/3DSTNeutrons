//File: NeutronCand.cpp
//Brief: NeutronCand is an Analyzer that plots quantities related to neutron candidates found as MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/NeutronCand.h"
#include "app/Factory.cpp"

//util includes
#include "ROOT/Base/TFileSentry.h"

namespace
{
  //TODO: If I am going to reuse this, put it into a class or at least a common header
  void Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::vector<int>& ids)
  {
    for(const auto& traj: trajs)
    {
      if(traj.ParentId == parent)
      {
        const int id = traj.TrackId;
        ids.push_back(id);
        Descendants(id, trajs, ids);
      }
    }
  }

  const TG4Trajectory& Matriarch(const TG4Trajectory& child, const std::vector<TG4Trajectory>& trajs)
  {
    if(child.ParentId == -1) return child;
    return Matriarch(trajs[child.ParentId], trajs);
  }
}

namespace ana
{
  NeutronCand::NeutronCand(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fClusters(*(config.Reader), "MergedClusters"), fMinEnergy(2.0)
  {
    fCandidateEnergy = config.File->make<TH1D>("CandidateEnergy", "Energy Specturm of Neutron Candidates;Energy [MeV];Events",
                                               150, 0, 1000);
    fCandPerNeutron = config.File->make<TH1D>("CandPerNeutron", "Number of Candidates per FS Neutron;Neutron Candidates;Neutrons",
                                              50, 0, 20);
    fNCand = config.File->make<TH1I>("NCand", "Number of Candidates per Event;Neutron Candidates;Events", 20, 0, 20);
    fFSNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of FS Neutrons that Produced Candidates;Energy [MeV];Events",
                                               200, 0, 3000);
    fCauseEnergyVsCandEnergy = config.File->make<TH2D>("CauseEnergyVsCandEnergy", "KE of FS Neutrons versus Energies of their Candidates;Candidate Energy [MeV];"
                                                                                  "FS Neutron KE [MeV];FS Neutrons", 150, 0, 1000, 200, 0, 3000);
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
          //fFSNeutronEnergy->Fill(trajs[part.TrackId].InitialMomentum.E()-trajs[part.TrackId].InitialMomentum.Mag());

          std::vector<int> descend;
          Descendants(part.TrackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          descend.push_back(part.TrackId);
          for(const auto& id: descend) TrackIDsToFS[id] = part.TrackId;
        }
      }
    }

    const double vtxBoxWidth = 10.;
    fNCand->Fill(std::count_if(fClusters.begin(), fClusters.end(), [&vtxBoxWidth, &TrackIDsToFS, &trajs](const auto& cand)
                                                                   {
                                                                     const auto& neutron = trajs[TrackIDsToFS[cand.TrackIDs.front()]];
                                                                     const auto vtxDiff = (cand.Position - neutron.Points[0].Position).Vect();
                                                                     return (vtxDiff.X() > vtxBoxWidth && 
                                                                             vtxDiff.Y() > vtxBoxWidth && 
                                                                             vtxDiff.Z() > vtxBoxWidth);
                                                                   }));

    std::map<int, std::vector<pers::MCCluster>> FSToCands; //Map of FS neutron to the candidates it is responsible for
    for(const auto& cand: fClusters)
    {
      //Don't count candidates near the vertex since my octree causes artifacts there.
      const auto vtxDiff = (cand.Position - trajs[TrackIDsToFS[cand.TrackIDs.front()]].Points[0].Position).Vect(); //TODO: This might get ugly 
                                                                                                            //      with multiple vertices.
      if(vtxDiff.X() > vtxBoxWidth && vtxDiff.Y() > vtxBoxWidth && vtxDiff.Z() > vtxBoxWidth) //10cm^3 vertex box
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
                                                                //neutrons because of Fermi motion.  
      }
    }

    //fCandPerNeutron->Fill(fClusters.GetSize()/FSToCands.size());

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
      fCandPerNeutronVsNeutronKE->Fill(trajs[FS.first].InitialMomentum.E()-trajs[FS.first].InitialMomentum.Mag(), FS.second.size());
    }
  }

  REGISTER_PLUGIN(NeutronCand, plgn::Analyzer);
}
