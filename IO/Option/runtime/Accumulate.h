//File: Accumulate.h
//Brief: An opt::Accumulate specifies how to handle options from the command line.  Accumulate concatenates the value 
//       it is given to the previous values with a space each time its key is specified.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_ACCUMULATE_H
#define OPT_ACCUMULATE_H

//util includes
#include "Policy.h"

//c++ includes that I can't get by without
#include <string>

namespace opt
{
  class Accumulate: public Policy
  {
    public:
      Accumulate(): Policy("<value>", "Accumulates all values passed to this option.  May be specified more than once.  ") {}
      virtual ~Accumulate() = default;

      std::string FoundFirst(const std::string& cmdLine) const override;
      std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const override;
      std::string FoundFirst() const override;
      std::string FoundAgain(const std::string& prev) const override;

  };
}

#endif //OPT_ACCUMULATE_H 
