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
}

namespace ana
{
  NeutronCand::NeutronCand(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fClusters(*(config.Reader), "MergedClusters"), fMinEnergy(2.0)
  {
    fCandidateEnergy = config.File->make<TH1D>("CandidateEnergy", "Energy Specturm of Neutron Candidates;Energy [MeV];Events",
                                               150, 0, 1000);
    fCandPerNeutron = config.File->make<TH1D>("CandPerNeutron", "Number of Candidates per FS Neutron;Neutron Candidates;Neutrons",
                                              50, 0, 5);
    fNCand = config.File->make<TH1I>("NCand", "Number of Candidates per Event;Neutron Candidates;Events", 5, 0, 5);
    fFSNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "Energies of FS Neutrons that Produced Candidates;Energy [MeV];Events",
                                               300, 0, 2000);
    fCauseEnergyVsCandEnergy = config.File->make<TH2D>("CauseEnergyVsCandEnergy", "Energies of FS Neutrons versus Energies of their Candidates;Candidate Energy [MeV];"
                                                                                  "FS Neutron Energy [MeV];FS Neutrons", 150, 0, 1000, 300, 0, 2000);
    fCandAngleWRTCause = config.File->make<TH1D>("CandAngleWRTCause", "Angle of Candidate w.r.t. InitialMomentum of FS Neutron;"
                                                                      "#Delta#theta_{Cand};Candidates", 200, 0, 3.1415926535897932384626433832);
    fDistFromVtx = config.File->make<TH1D>("DistFromVertex", "Distance of Closest Candidate to Vertex per FS Neutron;Distance [mm];FS Neutrons", 500, 0, 3000);
  }

  void NeutronCand::DoAnalyze()
  {
    std::map<int, int> TrackIDsToFS; //Map from TrackIDs to FS neutron
    const auto trajs = fEvent->Trajectories;

    std::cout << "Looping over " << fEvent->Primaries.size() << " primaries.\n";
    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        std::cout << "Looking at TrackID " << part.TrackId << " which is a " << part.Name << "\n";
        if(part.PDGCode == 2112 && part.Momentum.E() > fMinEnergy)
        {
          std::vector<int> descend;
          Descendants(part.TrackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          for(const auto& id: descend) TrackIDsToFS[id] = part.TrackId;
          std::cout << "Entered " << descend.size() << " descendants for TrackID " << part.TrackId << "\n";
        }
      }
    }

    fNCand->Fill(fClusters.GetSize());

    /*std::map<int, std::vector<pers::MCCluster>> FSToCands; //Map of FS neutron to the candidates it is responsible for
    for(const auto& cand: fClusters)
    {
      fCandidateEnergy->Fill(cand.Energy);
     
      std::set<int> FSIds; //The TrackIDs of FS neutrons responsible for this candidate.  Should almost always be only 1.
      for(const auto& id: cand.TrackIDs) FSIds.insert(TrackIDsToFS[id]); //TODO: The lookup of TrackIDsToFS is giving a default-constructed int here
      
      for(const int neutronID: FSIds) //For each FS neutron TrackID
      {
        FSToCands[neutronID].push_back(cand);
        const auto& neutron = trajs[neutronID];
        std::cout << "neutronID is " << neutronID << "\n";

        fCauseEnergyVsCandEnergy->Fill(cand.Energy, neutron.InitialMomentum.E());

        //Figure out the angle of the candidate w.r.t. the FS neutron's initial momentum
        const auto candVec = (cand.Position-neutron.Points[0].Position).Vect();

        fCandAngleWRTCause->Fill(std::acos(candVec.Unit().Dot(neutron.InitialMomentum.Vect().Unit())));
      }
    }

    fCandPerNeutron->Fill(fClusters.GetSize()/FSToCands.size());

    for(const auto& pair: FSToCands) fFSNeutronEnergy->Fill(trajs[pair.first].InitialMomentum.E());

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
    }*/
  }

  REGISTER_PLUGIN(NeutronCand, plgn::Analyzer);
}
