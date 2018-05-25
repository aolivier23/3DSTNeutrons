//File: Help.h
//Brief: Dedicated policy to support the "-h" and "--help" options.  Could be used for things like "--help-<other-option>" as well.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Policy.h"

namespace opt
{
  class Help: public Policy
  {
    public:
      Help(); 
      virtual ~Help() = default;
    
      virtual std::string FoundFirst(const std::string& cmdLine) const override;
      virtual std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const override;
      virtual std::string FoundFirst() const override;
      virtual std::string FoundAgain(const std::string& prev) const override;
  };
}
