//File: ExactlyOnce.cxx
//Brief: A command line option tht expects its key to be specified exactly once.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "ExactlyOnce.h"

namespace opt
{
  ExactlyOnce::ExactlyOnce(): Policy("<value>", "Must be specified exactly once.  ")
  {
  }

  std::string ExactlyOnce::FoundFirst(const std::string& cmdLine) const
  {
    return cmdLine;
  }

  std::string ExactlyOnce::FoundAgain(const std::string& cmdLine, const std::string& /*prev*/) const
  {
    throw exception("Argument found a second time, with value "+cmdLine+", that should only be found once"); 
  }

  std::string ExactlyOnce::FoundFirst() const
  {
    throw exception("Expected a value");
  }

  std::string ExactlyOnce::FoundAgain(const std::string& /*prev*/) const
  {
    throw exception("Expected a value");
  }
}
