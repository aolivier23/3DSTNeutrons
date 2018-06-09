//File: Accumulate.cxx
//Brief: An Accumulate Policy represents a command line option that can be specified multiple times.  Each new value specified 
//       is concatenated onto the previous values with a space.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Accumulate.h"

namespace opt
{
  std::string Accumulate::FoundFirst(const std::string& cmdLine) const
  {
    return cmdLine;
  }

  std::string Accumulate::FoundAgain(const std::string& cmdLine, const std::string& prev) const
  {
    return prev+" "+cmdLine;
  }

  //This Policy expects to be given a value always, so throw an exception if not given a value
  std::string Accumulate::FoundFirst() const
  {
    throw exception("Expected a value");
  }

  std::string Accumulate::FoundAgain(const std::string& ) const
  {
    throw exception("Expected a value");
  }
}
