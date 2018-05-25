//File: Exists.cxx
//Brief: A Exists option counts how many times it has been specified.  It expects to never get a value from the command line.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Exists.h"

namespace opt
{
  //This Policy expects to never be given a value
  std::string Exists::FoundFirst(const std::string& cmdLine) const
  {
    throw exception("Got value "+cmdLine+" where none was expected");
  }

  std::string Exists::FoundAgain(const std::string& cmdLine, const std::string& /*value*/) const
  {
    throw exception("Got value "+cmdLine+" where none was expected");
  }

  std::string Exists::FoundFirst() const
  {
    return "true";
  }

  std::string Exists::FoundAgain(const std::string& value) const
  {
    throw exception("Option was found twice that should only be specified once");
  }
} 
