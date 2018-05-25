//File: Exists.h
//Brief: An opt::Exists specifies how to handle options from the command line.  Exists represents an option 
//       that can be specified only once but should not be assigned to a value.  Throws if it gets a value.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_EXISTS_H
#define OPT_EXISTS_H

//util includes
#include "Policy.h"

//c++ includes that I can't get by without
#include <string>

namespace opt
{
  class Exists: public Policy
  {
    public:
      Exists(): Policy("", "Records whether this key was specified.  ") {}
      virtual ~Exists() = default;

      std::string FoundFirst(const std::string& cmdLine) const override;
      std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const override;
      std::string FoundFirst() const override;
      std::string FoundAgain(const std::string& prev) const override;
  };
} 

#endif //OPT_COUNTER_H
