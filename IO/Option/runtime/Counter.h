//File: Counter.h
//Brief: An opt::Counter specifies how to handle options from the command line.  Counter represents an option 
//       that can be specified multiple times but should not be assigned to a value.  Throws if it gets a value.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_COUNTER_H
#define OPT_COUNTER_H

//util includes
#include "Policy.h"

//c++ includes that I can't get by without
#include <string>

namespace opt
{
  class Counter: public Policy
  {
    public:
      Counter(): Policy("", "Counts how many times this key was specified.  ") {}
      virtual ~Counter() = default;

      std::string FoundFirst(const std::string& cmdLine) const override;
      std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const override;
      std::string FoundFirst() const override;
      std::string FoundAgain(const std::string& prev) const override;
  };
} 

#endif //OPT_COUNTER_H
