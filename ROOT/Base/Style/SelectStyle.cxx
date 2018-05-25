//File: SelectStyle.cxx
//Brief: Selects one of the style files in this directory to run based on a string.
//Author: Andrew Olivier aolivier@ur.rochester.edu
//TODO: I could probably combine CMake and/or preprocessor macros to have this file updated automatically.  

//local includes
#include "StandardStyle.cxx"
#include "DebugStyle.cxx"
//TODO: Include other useful styles here and below.

//c++ includes
#include <iostream> //For user feedback when using no style

namespace util
{
  void SelectStyle(const std::string& style)
  {
    if(style == "debug")
    {
      util::SetDebugStyle();
    }
    else if(style == "standard")
    {
      util::SetStandardStyle();
    }
    else
    {
      std::cerr << "Did not get a recognized style option, so not setting any style.\n";
    }
  }
}
