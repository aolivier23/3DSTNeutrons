//File: PrintCmdLine.cxx
//Brief: Logic to print out the command line to a file such that it can be made into an executable and a job be repeated.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_IO_OPTION_PRINTCMDLINE_CXX
#define UTIL_IO_OPTION_PRINTCMDLINE_CXX

//c++ includes
#include <fstream>
#include <string>

namespace opt
{
  void PrintCmdLine(int argc, char** argv, const std::string& producedFileName)
  {
    std::ofstream oFile("JobThatMade_"+producedFileName.substr(0, producedFileName.find_last_of("."))+".sh");
    oFile << "#!/bin/bash\n";
    for(int arg = 0; arg < argc; ++arg)
    {
      std::string argStr = argv[arg];
      if(argStr.find_first_of("(),'\"\\` #$&'") != std::string::npos)  //If this argument contains anything that Bash might consider 
                                                                       //a special character, enclose it in ''
      {
        argStr.insert(0, "'");
        argStr.push_back('\'');
      }
      oFile << argStr << " ";
    }
  }
}

#endif //UTIL_IO_OPTION_PRINTCMDLINE_CXX
