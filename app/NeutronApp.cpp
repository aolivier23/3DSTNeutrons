//File: NeutronApp.cpp
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
/*#include "IO/Option/runtime/CmdLine.h"
#include "IO/Option/runtime/Options.h"
#include "IO/Option/runtime/ExactlyOnce.h"
#include "IO/Option/runtime/Exists.h"
#include "IO/Option/runtime/Accumulate.h"
#include "IO/Option/PrintCmdLine.cxx"*/

//yaml-cpp includes
#include "yaml-cpp/yaml.h"

//ROOT includes
#include "TChain.h"
#include "TTreeReader.h" 
#include "TTreeReaderValue.h"
#include "TFile.h"
#include "TGeoManager.h" //For gGeoManager

//c++ includes
#include <iostream>

int main(int argc, const char** argv)
{
  try
  {
    //Read command line for application options.  Ignore keys that weren't requested for now since plugins might want them.  
    /*opt::CmdLine cmdLine("Runs reconstruction and analysis plugins over the entries in edepsim input files.  Plugins have the opportunity "
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

    const auto options = cmdLine.Parse(argc, argv);*/

    std::vector<std::string> inFiles;
    std::string configFiles; //Accumulate the content of all configuration files into this string

    //Parse the command line
    for(int pos = 1; pos < argc; ++pos) //Argument 0 is the path to this executable
    {
      const std::string arg(argv[pos]);
      if(arg.find(".yaml") != std::string::npos) 
      {
        std::ifstream input; //(arg);
        input.open(arg);
        if(!input) //If couldn't find configuration file in relative or absolute path, check in installation directory
        {
          std::cout << "Failed to read file " << arg << ", so looking in package installation directory...\n";
          input.clear();
          std::string path = "";
          const auto value = getenv("THREEDSTNEUTRONS_CONF_PATH");
          if(value != nullptr) 
          {
            path = value;
            path += "yaml/";
          }
          input.open(path+arg);
        }
        if(!input) //Assume that file doesn't exist
        {
          std::cerr << "Failed to find configuration file named " << arg << "\n";
          return 8;
        }

        configFiles.append(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
      }
      else if(arg.find(".root") != std::string::npos) inFiles.push_back(arg);
      else 
      {
        std::cerr << "Got command line argument that is neither a ROOT file nor a configuration file:" << arg << "\n";
        return 7;
      }
      //TODO: Regular expression support via RegexFiles?
    }

    YAML::Node config = YAML::Load(configFiles);
    //std::cout << "Configuration state is:\n" << config << "\n";

    //TODO: Include directives?  For now, specify "include"d files on the command line before other files need them.  

    //Look for options for the application first
    long int nEvents = -1; //placeholder value
    if(config["app"])
    {
      const auto& appOpt = config["app"];
      if(appOpt["source"])
      {
        const auto& source = appOpt["source"];
        const auto& files = source["files"];
        if(files) for(auto name = files.begin(); name != files.end(); ++name) inFiles.push_back(name->as<std::string>());
        //TODO: Support for regular expressions via RegexFiles?

        if(source["NEvents"]) nEvents = source["NEvents"].as<long int>();
      }
    }

    //Validate configuration so far and prepare to read files
    if(inFiles.empty())
    {
      std::cerr << "No input files found, so not doing anything.\n";
      return 6;
    }

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

    TFile* outFile = nullptr;
    TTree* outTree = nullptr;

    //Find all algorithms from the configuration document.  
    std::vector<std::unique_ptr<plgn::Reconstructor>> recoAlgs;

    if(config["reco"])
    {
      //Only create an output file if there are Reconstructors being run
      //Create a copy of the structure of the input tree
      outFile = TFile::Open(config["reco"]["OutputName"].as<std::string>().c_str(), "CREATE"); 
      if(!outFile)
      {
        std::cerr << "Could not create a new file called " << outFile->GetName() << " to write out reconstructed events.\n";
        return 3;
      }

      outTree = inTree->CloneTree(0);
      outTree->SetDirectory(outFile);

      plgn::Reconstructor::Config recoConfig;
      recoConfig.Input = &inReader;
      recoConfig.Output = outTree;

      const auto& recos = config["reco"]["algs"];
      auto& recoFactory = plgn::Factory<plgn::Reconstructor>::instance();
      for(auto reco =  recos.begin(); reco != recos.end(); ++reco)
      {
        recoConfig.Options = reco->second;
        auto recoAlg = recoFactory.Get(reco->first.as<std::string>(), recoConfig);
        if(recoAlg)
        {
          recoAlgs.push_back(std::move(recoAlg));
        }
        else std::cerr << "Could not find Reconstructor algorithm " << reco->first << "\n";
      }
    }
    else std::cout << "No Reconstructors specified, so not creating an output file.\n";
 
    //Parameters from the command line
    //const auto inFiles = util::RegexFilesPath<std::string>(options["--regex"], options["--path"]);

    //Get analysis plugins
    std::vector<std::pair<std::string, std::unique_ptr<plgn::Analyzer>>> anaAlgs;    

    std::unique_ptr<util::TFileSentry> anaFile(nullptr); 
    if(config["analysis"]) 
    { 
      anaFile.reset(new util::TFileSentry(config["analysis"]["FileName"].as<std::string>().c_str())); 
      //Only create histogram file if I am running analysis plugins                                                              
      util::SelectStyle(config["analysis"]["style"].as<std::string>()); 
      plgn::Analyzer::Config anaConfig;
      anaConfig.File = anaFile.get();
      anaConfig.Reader = &inReader;
      //anaConfig.Options = &options;
  
      const auto& anas = config["analysis"]["algs"];
      auto& anaFactory = plgn::Factory<plgn::Analyzer>::instance();
      for(auto ana = anas.begin(); ana != anas.end(); ++ana)
      {
        anaConfig.Options = ana->second;
        anaFile->cd(ana->first.as<std::string>());
        auto anaAlg = anaFactory.Get(ana->first.as<std::string>(), anaConfig);
        if(anaAlg) anaAlgs.emplace_back(ana->first.as<std::string>(), std::move(anaAlg));
        else std::cerr << "Could not find Analyzer algorithm " << ana->first << "\n";
      }
    }
    else std::cout << "No Analyzers specified, so not creating a histogram file.\n";

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

      if(outTree) inTree->CopyAddresses(outTree);
      else std::cout << "There is no output tree, so not copying addresses.\n"; //TODO: This is only debugging output.  Remove it from release builds?
      inReader.SetTree(inTree);

      std::cout << "Processing file " << file << "\n";

      for(const auto entry: inReader)
      {
        if(entry == nEvents) break;

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
          if(!outTree) std::cerr << "Did some reconstruction, but output TTree has not been created!\n"; //TODO: This is only debugging output.  Remove it from release builds?
          outTree->Fill();
        }
        //else std::cout << "No reconstruction objects to save for event " << entry << "\n";

        //Next, call analysis plugins
        for(const auto& ana: anaAlgs) 
        {
          //TODO: Change to directory for this analyzer in case make is called during Analyze.  This might be an indication that I need to rethink
          //      TFileSentry.
          anaFile->cd(ana.first);
          ana.second->Analyze();
        }

        if(entry%100 == 0 || entry < 100) std::cout << "Finished processing event " << entry << "\n";

        //TODO: Use gGeoManager in plugins for now, but consider retrieving TGeoManager from current file instead.  
      }
    }
   
    //Write out the reconstruced TTree if there was any reconstruction done.  
    if(outFile)
    {
      if(!outTree) std::cerr << "Output file was created, but there is no output TTree!\n";
      outFile->cd();
      outTree->Write(); //TODO: Create a friend TTree instead of copying the entire file to reduce size of reconstruction.  See Guang's port of 
                        //      T2K FGD code for example of friend trees.

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
    else std::cout << "No output file created, so nothing to write.  Histograms written when main() ends.\n";
  }
  catch(const std::exception& e)
  {
    std::cerr << "Caught STL exception:\n" << e.what() << "\n";
    return 4;
  }
  catch(const util::exception& e)
  {
    std::cerr << e.what() << "\n";
    return 5;
  }
}
