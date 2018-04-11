//File: NeutronTOF.cpp
//Brief: NeutronTOF is an Analyzer that plots quantities related to neutron candidates found as MCHits.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EDepNeutrons includes
#include "ana/NeutronTOF.h"
#include "app/Factory.cpp"
#include "alg/TruthFunc.h"

//util includes
#include "ROOT/Base/TFileSentry.h"
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"

//c++11 includes
#include <chrono>

namespace plgn
{
  //Register command line options
  template <>
  void RegCmdLine<ana::NeutronTOF>(opt::CmdLine& opts)
  {
    opts.AddKey("--hit-alg", "Name of the branch of pers::MClusters to analyze.  Usually the name of the cluster-making algorithm.", "MergedClusters");
    opts.AddKey("--time-res", "Toy timing resolution of a 3DST in ns.  Used for binning and smearing hit times in the NeutronTOF algorithm.", "0.7");
  }
}

namespace ana
{
  NeutronTOF::NeutronTOF(const plgn::Analyzer::Config& config): plgn::Analyzer(config), fHits(*(config.Reader), (*(config.Options))["--hit-alg"].c_str()), 
                                                                fGen(std::chrono::system_clock::now().time_since_epoch().count()), 
                                                                fGaus(0., config.Options->Get<double>("--time-res"))
  {
    const float timeMax = 100., distMax = 5000.;
    const size_t nTimeBins = timeMax/config.Options->Get<double>("--time-res"), nDistBins = distMax/10.0; //1.0cm bins
    fNeutronHitTime = config.File->make<TH1D>("NeutronHitTime", "Time of First Hit from a FS Neutron;Time [ns];Visible FS Neutrons", nTimeBins, 0, timeMax); 

    fNeutronTimeVersusDist = config.File->make<TH2D>("NeutronTimeVersusDist", "Time of First Hit from a FS Neutron Versus Distance;Distance [mm];Time [ns]",
                                                     nDistBins, 0, distMax, nTimeBins, 0, timeMax);                                                     

    fNeutronTOFEnergy = config.File->make<TH1D>("NeutronTOFEnergy", "Kinetic Energy from TOF and Distance to First Hit for FS Neutrons;Energy [MeV]", 300, 0., 1000.);

    fNeutronEResidual = config.File->make<TH1D>("NeutronEResidual", "Relative Error in Neutron Energy from TOF;#frac{E_{TOF}-E_{True}}{E_{True}}", 300, -1., 1.);

    fBeta = config.File->make<TH1D>("Beta", "Velocity Ratio for Closest Hit to Each FS Neutron;#frac{v}{c}", 50, 0, 1.);
  }

  void NeutronTOF::DoAnalyze()
  {
    std::map<int, int> TrackIDsToFS; //Map from TrackIDs to FS neutron
    const auto trajs = fEvent->Trajectories;

    for(const auto& vertex: fEvent->Primaries)
    {
      for(const auto& part: vertex.Particles)
      {
        if(part.PDGCode == 2112) //&& part.Momentum.E() - part.Momentum.Mag() > fMinEnergy)
        {
          std::set<int> descend;
          truth::Descendants(part.TrackId, trajs, descend); //Fill descend with the TrackIDs of part's descendants
          descend.insert(part.TrackId);
          for(const auto& id: descend) TrackIDsToFS[id] = part.TrackId; 
        }
      }
    }

    for(const auto& vert: fEvent->Primaries)
    {
      auto vertPos = vert.Position;
      for(const auto& part: vert.Particles)
      {
        if(part.Name == "neutron")
        {
          //Find the closest hit to vert for this neutron
          auto closest = fHits.end();
          for(auto iter = fHits.begin(); iter != fHits.end(); ++iter)
          {
            auto& hit = *iter;
            if(std::find_if(hit.TrackIDs.begin(), hit.TrackIDs.end(), 
                            [&TrackIDsToFS, &part](const int id) { return TrackIDsToFS[id] == part.TrackId; }) 
               != hit.TrackIDs.end())
            {
              if(closest == fHits.end()) closest = iter;
              else if((hit.Position - vertPos).Vect().Mag() < ((*closest).Position - vertPos).Vect().Mag()) closest = iter;
            }
          }

          //If I found a first hit for this FS neutron
          if(closest != fHits.end())
          {
            const auto diff = ((*closest).Position - vert.Position); 
            const auto smear = fGaus(fGen);
            const double deltaT = diff.T() - (vert.Position.T()+smear); //Smear vertex time values since I'm using the true vertex for now
            const double dist = diff.Vect().Mag(); //Distance is already smeared by virtue of the geometry I am using for hit-making
            fNeutronHitTime->Fill(deltaT);
            fNeutronTimeVersusDist->Fill(dist, deltaT);

            if(deltaT > 3.0 && dist > 10.) //in ns and mm respectively
            {
              const float c = 299.792; //Speed of light = 300 mm/ns
              const auto beta = dist/deltaT/c; 
              const auto mass = 939.56563; //MeV/c^2
              const auto energy = mass/std::sqrt(1.-beta*beta); //E = gamma * mc^2 
              fBeta->Fill(beta);
              fNeutronTOFEnergy->Fill(energy-mass); //TODO: Error bars

              const auto trueE = part.Momentum.E();
              fNeutronEResidual->Fill((energy - trueE)/trueE);

              if(beta < 0.02) 
              {
                std::cout << "Beta is < 0.02: " << beta << ".  Distance is " << dist
                          << "\nTime difference is " << deltaT << "\n"
                          << "Interaction time is " << vert.Position.T() << "\n"
                          << "Closest hit time is " << (*closest).Position.T() << "\n"
                          << "Smeared vertex time by " << smear << "\n"
                          << "closest->Position is (" << (*closest).Position.X() << ", " << (*closest).Position.Y() << ", " 
                          << (*closest).Position.Z() << ")\n"
                          << "Vertex is (" << vert.Position.X() << ", " << vert.Position.Y() << ", " << vert.Position.Z() << ")\n"
                          << "EventID is " << fEvent->EventId << "\n";
              }

              //const auto uncert = beta*m/gamma/gamma/gamma/deltaT*std::sqrt(1.*1./c/c+0.7*0.7*beta*beta);
            } //If time difference and distance are above minimal resolution
          } //If there is a visible hit for this FS neutron
        } //If this particle is a neutron
      } //For each particle in this primary vertex
    } //For each primary vertex
  }

  REGISTER_PLUGIN(NeutronTOF, plgn::Analyzer);
}
