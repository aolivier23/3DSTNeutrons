//File: ExactlyOnce.h
//Brief: An opt::ExactlyOnce specifies how to handle options from the command line.  ExactlyOnce looks for any 
//       string that is only set once.  Throws on trying to set option value again.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_EXACTLYONCE_H
#define OPT_EXACTLYONCE_H

//util includes
#include "Policy.h"

//c++ includes that I can't get by without
#include <string>

namespace opt
{
  class ExactlyOnce: public Policy
  {
    public:
      ExactlyOnce();
      virtual ~ExactlyOnce() = default;

      virtual std::string FoundFirst(const std::string& cmdLine) const override;
      virtual std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const override;
      virtual std::string FoundFirst() const override;
      virtual std::string FoundAgain(const std::string& prev) const override;

  };
} 

#endif //OPT_EXACTLYONCE_H
