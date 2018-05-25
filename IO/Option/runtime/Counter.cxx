//File: Counter.cxx
//Brief: A Counter option counts how many times it has been specified.  It expects to never get a value from the command line.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Counter.h"

namespace opt
{
  //This Policy expects to never be given a value
  std::string Counter::FoundFirst(const std::string& cmdLine) const
  {
    throw exception("Got value "+cmdLine+" where none was expected");
  }

  std::string Counter::FoundAgain(const std::string& cmdLine, const std::string& /*value*/) const
  {
    throw exception("Got value "+cmdLine+" where none was expected");
  }

  std::string Counter::FoundFirst() const
  {
    return "1";
  }

  std::string Counter::FoundAgain(const std::string& value) const
  {
    size_t val = atoi(value.c_str());
    return std::to_string(val+1);
  }
} 
