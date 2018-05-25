//File: RuntimeOptExample.cpp
//Brief: Example of how to the the CmdLine runtime option parser to get information from the command line.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//opt includes
#include "CmdLine.h" 
#include "ExactlyOnce.h"
#include "Counter.h"
#include "Accumulate.h"
#include "Exists.h"
#include "Options.h"

//c++ includes
#include <iostream>
#include <vector>
//#include <fstream>

int main(int argc, char** argv)
{
  try
  {
    //Set up the command line parser
    opt::CmdLine cmdLine("Demonstrates how to use the opt::CmdLine class to parse the command line at runtime.  Accesses various options.");
    //                         key      help text                                                                        default value
    cmdLine.AddKey(            "--print", "Prints <text> to STDOUT.  ",                                     "Hello World!");
    cmdLine.AddKey<opt::Counter>(   "-v",     "Specify verbosity.  Can be given up to three times.  ",                    "0");
    cmdLine.AddKey<opt::Accumulate>("--file",  "Specify one or more files.  Can be called multiple times.  ");
    cmdLine.AddKey<opt::Exists>("--exists", "Check whether this option was specified.  ", "false");
    const auto options = cmdLine.Parse(argc, argv);

    const size_t verbosity = options.Get<int>("-v");
    if(verbosity > 0) std::cout << "option -v was specified " << verbosity << " times.\n";
    if(verbosity > 2) std::cout << "About to print value for option --print:\n";
    std::cout << options["--print"] << "\n";

    const bool exists = options.Get<bool>("--exists");
    if(exists) std::cout << "Option --exists was specified.\n";

    if(verbosity > 1) std::cout << "File names are:\n";
    const auto files = options.Get<std::vector<std::string>>("--file");
    for(const auto& file: files) std::cout << file << "\n";
    
    //Print the command line to rerun this job to a bash script
    //TODO: This doesn't work with Counter policies
    //std::ofstream commands("commands.sh");
    //commands << options.PrintToBash();
  }
  catch(const opt::CmdLine::exception& e)
  {
    std::cerr << e.what() << "\n";
    return -1;
  }
  catch(const opt::Options::exception& e)
  {
    std::cerr << e.what() << "\n";
    return -2;
  }

  return 0;
}
