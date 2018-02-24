//File: EdepApp.cpp
//Brief: An EdepApp creates user-specified reconstruction and analysis plugins and processes a user-specified edepsim output file.  
//       The reconstruction plugins are given read-only access to each entry from the old TTree and write access to 
//       an entry in the new TTree. Reconstruction objects are written directly to the edepsim tree, as opposed to 
//       included in the TG4Event.  
//
//       Analysis plugins are given read-only access to the read in TTree after the reconstruction plugins have processed it.  
//       Since the output tree is a clone of the input tree, analysis plugins will see any changes the reconstruction plugins 
//       made.   
//Author: Andrew Olivier aolivier@ur.rochester.edu

//edepsim includes
#include "TG4Event.h"

//Plugin includes
#include "app/Factory.cpp"
#include "ana/Analyzer.h"
#include "reco/Reconstructor.h"

//util includes
#include "IO/File/RegexFiles.cxx"
#include "ROOT/Base/Style/SelectStyle.cxx"
#include "ROOT/Base/TFileSentry.h"

//Option parser includes
#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"
#include "IO/Option/runtime/Exists.h"
#include "IO/Option/runtime/Accumulate.h"
#include "IO/Option/PrintCmdLine.cxx"

//ROOT includes
#include "TChain.h"
#include "TTreeReader.h" 
#include "TTreeReaderValue.h"
#include "TFile.h"
#include "TGeoManager.h" //For gGeoManager

//c++ includes
#include <iostream>

int main(int argc, char** argv)
{
  try
  {
    //Read command line for application options.  Ignore keys that weren't requested for now since plugins might want them.  
    opt::CmdLine cmdLine("Runs reconstruction and analysis plugins over the entries in edepsim input files.  Plugins have the opportunity "
                         "to request additional objects not in normal edepsim TTrees.");
    cmdLine.AddKey("--regex", "Regular expression for files to read.");
    cmdLine.AddKey("--path", "Path to files to be read.");
    cmdLine.AddKey<opt::Accumulate>("--reco", "List of reconstruction plugins to run.", "");
    cmdLine.AddKey<opt::Accumulate>("--ana", "List of analysis plugins to run.", "");
    cmdLine.AddKey("--style", "Specify how to configure the style for drawing ROOT objects.  Options are: none, debug, standard.", "debug");
    cmdLine.AddKey("--reco-file", "The name of the file to which the reconstructed TTree will be saved.");
    cmdLine.AddKey("--n-events", "The maximum number of events to process.", "-1");

    //Get command line options specified by plugins at compile-time.
    auto& recoFactory = plgn::Factory<plgn::Reconstructor>::instance();
    recoFactory.RegCmdLine(cmdLine);
    auto& anaFactory = plgn::Factory<plgn::Analyzer>::instance();
    anaFactory.RegCmdLine(cmdLine);

    const auto options = cmdLine.Parse(argc, argv);

    //Print out the command line so this job can be repeated easily
    opt::PrintCmdLine(argc, argv, options["--reco-file"]);

    //Parameters from the command line
    const auto inFiles = util::RegexFilesPath<std::string>(options["--regex"], options["--path"]);
    const std::string outFileName;

    //Set up to read files
    auto inFile = TFile::Open(inFiles.begin()->c_str(), "READ");
    if(!inFile)
    {
      std::cerr << "Failed to open file " << inFiles.begin()->c_str() << " for reading, so quitting.\n";
      return 1;
    }

    auto inTree = (TTree*)inFile->Get("EDepSimEvents");
    if(!inTree)
    {
      std::cerr << "File " << inFile->GetName() << " did not have a TTree named EDepSimEvents, so it is not an edepsim input file.\n";
      return 2;
    }         

    TTreeReader inReader(inTree);
    TTreeReaderValue<TG4Event> event(inReader, "Event");  

    //Create a copy of the structure of the input tree
    auto outFile = TFile::Open(options["--reco-file"].c_str(), "CREATE");
    if(!outFile)
    { 
      std::cerr << "Could not create a new file called " << options["--reco-file"] << " to write out reconstructed events.\n";
      return 3;
    }

    auto outTree = inTree->CloneTree(0); 
    outTree->SetDirectory(outFile);
    
    //Get reconstruction plugins
    plgn::Reconstructor::Config config;
    config.Input = &inReader;
    config.Output = outTree;
    config.Options = &options;

    const auto recos = options.Get<std::vector<std::string>>("--reco");
    std::vector<std::unique_ptr<plgn::Reconstructor>> recoAlgs;
    for(const auto& reco: recos)
    {
      if(reco == "") continue; //TODO: Make sure this never makes it into an Accumulate option
      auto recoAlg = recoFactory.Get(reco, config);
      if(recoAlg) recoAlgs.push_back(std::move(recoAlg));
      else std::cerr << "Could not find Reconstructor algorithm " << reco << "\n";
    }
                                                                        
    //Get analysis plugins
    //TODO: Directory for each plugin?
    util::SelectStyle(options["--style"]);
    util::TFileSentry anaFile("histos_"+options["--reco-file"]);
    plgn::Analyzer::Config anaConfig;
    anaConfig.File = &anaFile;
    anaConfig.Reader = &inReader;
    anaConfig.Options = &options;

    const auto anas = options.Get<std::vector<std::string>>("--ana");
    std::vector<std::unique_ptr<plgn::Analyzer>> anaAlgs;
    for(const auto& ana: anas)
    {
      if(ana == "") continue; //TODO: Make sure this never makes it into an Accumulate option
      auto anaAlg = anaFactory.Get(ana, anaConfig);
      if(anaAlg) anaAlgs.push_back(std::move(anaAlg));
      else std::cerr << "Could not find Analyzer algorithm " << ana << "\n";
    }

    for(const auto& file: inFiles)
    {
      if(inFile) delete inFile; //Make sure previous file is closed.  
      inFile = TFile::Open(file.c_str(), "READ");
      if(!inFile) 
      {
        std::cerr << "Could not open file " << file << " for reading, so skipping it.\n";
        continue;
      }

      inTree = (TTree*)inFile->Get("EDepSimEvents");
      if(!inTree)
      {
        std::cerr << "Could not find TTree named EDepSimEvents in " << file << ", so skipping this file.\n";
        continue;
      } 

      gGeoManager = (TGeoManager*)inFile->Get("EDepSimGeometry");
      if(!gGeoManager)
      {
        std::cerr << "Could not find a TGeoManager named EDepSimGeometry in " << file << ", so skipping this file.\n";
        continue;
      }

      inTree->CopyAddresses(outTree);
      inReader.SetTree(inTree);

      std::cout << "Processing file " << file << "\n";

      for(const auto entry: inReader)
      {
        //First, call Reconstructor plugins
        bool foundReco = false;
        for(const auto& reco: recoAlgs)
        {
          foundReco = (reco->Reconstruct())?true:foundReco;
        }

        //If something was reconstructed, write to the output tree
        if(foundReco) 
        {
          inTree->GetEntry(entry); //TODO: Why does this work when SetBranchStatus() doesn't?  mysteriesOfTheUniverse.push_back(this)
          outTree->Fill();
        }

        //Next, call analysis plugins
        for(const auto& ana: anaAlgs) ana->Analyze();

        if(entry%100 == 0 || entry < 100) std::cout << "Finished processing event " << entry << "\n";
        if(entry > options.Get<size_t>("--n-events")) break; 

        //TODO: Use gGeoManager in plugins for now, but consider retrieving TGeoManager from current file instead.  
      }
    }
   
    //Write out the reconstruced TTree.  
    outFile->cd();
    outTree->Write();

    //Copy geometry and edepsim PassThru information from last file (?)
    //TODO: Copy from all files
    //For now, assuming that the geometry is the same in each file and the pass-thru information is an empty directory.  
    auto man = (TGeoManager*)inFile->Get("EDepSimGeometry");
    //auto passThru = (TDirectoryFile*)inFile->Get("DetSimPassThru"); //This has always been empty so far, so not copying it for now.
    outFile->cd();
    man->Write();
    //passThru->Write();
    
    outFile->Write(); //TODO: Is this necessary?
  }
  catch(const std::exception& e)
  {
    std::cerr << "Caught STL exception:\n" << e.what() << "\n";
    return 1;
  }
  catch(const util::exception& e)
  {
    std::cerr << e.what() << "\n";
  }
}
