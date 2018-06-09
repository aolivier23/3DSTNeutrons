//File: CandTOF.cpp
//Brief: CandTOF is an Analyzer that plots quantities related to neutron candidates found as MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/CandTOF.h"
#include "app/Factory.cpp"
#include "alg/TruthFunc.h"

//util includes
#include "ROOT/Base/TFileSentry.h"

//c++11 includes
#include <chrono>

/*namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::CandTOF>(opt::CmdLine& opts)
  {
    opts.AddKey("--cand-alg", "Name of the branch of pers::NeutronCands to analyze.  Usually the name of the cluster-making algorithm.", "CandFromPDF");
    opts.AddKey("--clust-alg", "Name of the branch of pers::MClusters to analyze.  Usually the name of the cluster-making algorithm.", "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the CandTOF algorithm.", "0.7");
  }
}*/

namespace ana
{
  CandTOF::CandTOF(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fCands(*(config.Reader), 
                                                                                  config.Options["CandAlg"].as<std::string>().c_str()), 
                                                                fClusters(*(config.Reader), 
                                                                          config.Options["ClusterAlg"].as<std::string>().c_str()),
                                                                fGen(std::chrono::system_clock::now().time_since_epoch().count()), 
                                                                fGaus(0., config.Options["TimeRes"].as<double>()), fPosRes(10.), 
                                                                fTimeRes(config.Options["TimeRes"].as<double>())
  {
    const float timeMax = 100., distMax = 5000.;
    const size_t nTimeBins = timeMax/config.Options["TimeRes"].as<double>(), nDistBins = distMax/10.0; //1.0cm bins
    fNeutronHitTime = config.File->make<TH1D>("NeutronHitTime", "Time of First Hit from a FS Neutron;Time [ns];Visible FS Neutrons", nTimeBins, 0, timeMax); 

    fNeutronTimeVersusDist = config.File->make<TH2D>("NeutronTimeVersusDist", "Time of First Hit from a FS Neutron Versus Distance;Distance [mm];Time [ns]",
                                                     nDistBins, 0, distMax, nTimeBins, 0, timeMax);                                                     

    fCandTOFEnergy = config.File->make<TH1D>("CandTOFEnergy", "Kinetic Energy from TOF and Distance to First Hit for FS Neutrons;Energy [MeV]", 300, 0., 1000.);

    fNeutronEResidual = config.File->make<TH1D>("NeutronEResidual", "Relative Error in Neutron Energy from TOF;#frac{E_{TOF}-E_{True}}{E_{True}}", 300, -1., 1.);

    fBeta = config.File->make<TH1D>("Beta", "Velocity Ratio for Closest Hit to Each FS Neutron;#frac{v}{c}", 50, 0, 1.);
   
    fTrueBeta = config.File->make<TH1D>("TrueBeta", "Initial Velocity Ratios for Visible FS Neutrons;#frac{v}{c}", 50, 0, 1.);

    fBetaRes = config.File->make<TH1D>("BetaRes", "How Different is Neutron Speed from c in #sigma_{#beta}s;#sigma_{#beta}s", 20, 0, 20);

    fFSNeutronEnergy = config.File->make<TH1D>("FSNeutronEnergy", "KE of FS Neutrons that Produced Candidates;Energy [MeV];Events",
                                               200, 0, 3000);

    fTotalEResidual = config.File->make<TH1D>("TotalEResidual", "How well does Total Energy from TOF Represent Total Neutron Energy?;"
                                                                "Neutron Energy [MeV];Events", 300., -1., 1.);

    //TODO: New plot(s)
  }

  void CandTOF::DoAnalyze()
  {
    std::map<int, int> TrackIDsToFS; //Map from TrackIDs to FS neutron
    const auto trajs = fEvent->Trajectories;

    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const int pdg = part.GetPDGCode();
        const int trackId = part.GetTrackId();
        #else
        const int pdg = part.PDGCode;
        const int trackId = part.TrackId;
        #endif

        if(pdg == 2112) //&& part.Momentum.E() - part.Momentum.Mag() > fMinEnergy)
        {
          std::set<int> descend;
          truth::Descendants(trackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          descend.insert(trackId);
          for(const auto& id: descend) TrackIDsToFS[id] = trackId; 
        }
      }
    }

    //for(const auto& vert: fEvent->Primaries)
    const auto& vert = fEvent->Primaries.front(); //TODO: Associate NeutronCands with vertices?
    {
      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
      const auto& vertPos = vert.GetPosition();
      #else
      const auto& vertPos = vert.Position;
      #endif

      //Require that this is a CC interaction.  I want to avoid the problem of NC vertex reconstruction for now since it may not be 
      //possible in many cases. 
      if(std::find_if(vert.Particles.begin(), vert.Particles.end(), [](const auto& part) 
                                                                    { 
                                                                      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                                                      const int pdg = part.GetPDGCode();
                                                                      #else
                                                                      const int pdg = part.PDGCode;
                                                                      #endif

                                                                      return pdg == 13 || pdg == -13 || pdg == 12 || pdg == -13; 
                                                                    })
         != vert.Particles.end())
      {
  
        //Since I'm not reconstructing neutrino vertices yet, smear the vertex time by timing resolution
        const auto smear = fGaus(fGen);
  
        double totalTOFE = 0.;
  
        for(const auto& cand: fCands)
        {
          const auto diff = cand.Start - vertPos;
          const double deltaT = diff.T() - smear; //Smear vertex time values since I'm using the true vertex for now
          const double dist = diff.Vect().Mag(); //Distance is already smeared by virtue of the geometry I am using for hit-making
          fNeutronHitTime->Fill(deltaT);
          fNeutronTimeVersusDist->Fill(dist, deltaT);
  
          //TODO: Tune deltaT cut?
          if(deltaT > fTimeRes && dist > fPosRes) //in ns and mm respectively
          {
            //Find the highest-energy FS neutron that contributed to the NeutronCand
            auto highestID = cand.TrackIDs.end();
            double highestE = 0.;
            for(auto trackID = cand.TrackIDs.begin(); trackID != cand.TrackIDs.end(); ++trackID)
            {
              const auto& traj = trajs[*trackID];

              #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
              if(strcmp(traj.GetName(), "neutron") && traj.GetInitialMomentum().E() > highestE) highestID = trackID;
              #else
              if(traj.Name == "neutron" && traj.InitialMomentum.E() > highestE) highestID = trackID;
              #endif
            }

            const auto& part = *(std::find_if(vert.Particles.begin(), vert.Particles.end(), [&highestID](const auto& part) 
                                                                                            {
                                                                                              #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS 
                                                                                              return part.GetTrackId() == *highestID; 
                                                                                              #else
                                                                                              return part.TrackId == *highestID;
                                                                                              #endif
                                                                                            }));

            const float c = 299.792; //Speed of light = 300 mm/ns
            const auto beta = dist/deltaT/c; 
            const auto mass = 939.56563; //MeV/c^2
            const auto energy = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2 
            totalTOFE += energy - mass;
            fBeta->Fill(beta);
            const double distUncert = fPosRes/dist; //relative uncertainty in distance for this energy "measurement"
            const double timeUncert = fTimeRes/deltaT; //relative uncertainty in time for this energy "measurement"
            const double betaUncert = beta*std::sqrt(distUncert*distUncert+timeUncert*timeUncert); //Uncertainty in beta
            fBetaRes->Fill((1-beta)/betaUncert); //By how many sigmas is beta different from 1?

            #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
            const double trueGamma = part.GetMomentum().E()/part.GetMomentum().Mag();
            #else
            const double trueGamma = part.Momentum.E()/part.Momentum.Mag();
            #endif

            fTrueBeta->Fill(std::sqrt(1.-1./trueGamma/trueGamma));
            fCandTOFEnergy->Fill(energy-mass); //TODO: Error bars

            #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
            const auto trueE = part.GetMomentum().E();
            #else
            const auto trueE = part.Momentum.E();
            #endif

            fNeutronEResidual->Fill((energy - trueE)/trueE);
  
            if(beta < 0.02)
            {
              std::cout << "Beta is < 0.02: " << beta << ".  Distance is " << dist
                        << "\nTime difference is " << deltaT << "\n"
                        << "Interaction time is " << vertPos.T() << "\n"
                        << "Closest hit time is " << cand.Start.T() << "\n"
                        << "Smeared vertex time by " << smear << "\n"
                        << "closest->Position is (" << cand.Start.X() << ", " << cand.Start.Y() << ", " 
                        << cand.Start.Z() << ")\n"
                        << "Vertex is (" << vertPos.X() << ", " << vertPos.Y() << ", " << vertPos.Z() << ")\n"
                        << "EventID is " << fEvent->EventId << "\n";
            }
  
            //const auto uncert = beta*m/gamma/gamma/gamma/deltaT*std::sqrt(1.*1./c/c+0.7*0.7*beta*beta);
          } //If time difference and distance are above minimal resolution
        } //For each NeutronCand in this event
        const double totalTrueE = std::accumulate(vert.Particles.begin(), vert.Particles.end(), 0., [](const auto sum, const auto& part)
                                                                                                    { 
                                                                                                      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
                                                                                                      if(!strcmp(part.GetName(), "neutron")) 
                                                                                                      return sum;
                                                                                                      return sum + part.GetMomentum().E() - 
                                                                                                             part.GetMomentum().Mag();
                                                                                                      #else
                                                                                                      if(part.Name != "neutron") return sum;
                                                                                                      return sum + part.Momentum.E() -
                                                                                                             part.Momentum.Mag(); 
                                                                                                      #endif
                                                                                                    });
  
        fTotalEResidual->Fill((totalTOFE-totalTrueE)/totalTrueE);
      } //If this vertex has a muon or an electron
    } //For each primary vertex
  }

  REGISTER_PLUGIN(CandTOF, plgn::Analyzer)
}
