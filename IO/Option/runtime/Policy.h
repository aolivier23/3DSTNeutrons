//File: Policy.h
//Brief: An opt::Policy specifies how to handle options from the command line.  In CmdLine, each key is paired with an opt::Policy. 
//       When that key is found on the command line, std::string opt::Policy::FoundFirst(const std::string& cmdLine) the first 
//       time that key is found.  If there is already a value for the key mapped to a Policy, 
//       std::string opt::Policy::FoundAgain(const std::string& cmdLine, const std::string& prev) will be called.  
//       Note that default values for keys are applied after parsing is completed.  Policy objects do not know about default values.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_POLICY_H
#define OPT_POLICY_H

//util includes
#include "Base/exception.h"

//c++ includes that I can't get by without
#include <string>

namespace opt
{
  class Policy
  {
    public:
      Policy(const std::string& usage, const std::string& details); //Define this in Policy.cxx to keep CMake happy while still 
                                                                    //needing to include only a .h file in other files.  
      virtual ~Policy() = default;

      virtual std::string FoundFirst(const std::string& cmdLine) const = 0; //For keys that expect a value
      virtual std::string FoundFirst() const = 0; //For keys that do not expect a value
      virtual std::string FoundAgain(const std::string& cmdLine, const std::string& prev) const = 0; //For keys that expect a value
      virtual std::string FoundAgain(const std::string& prev) const = 0; //For keys that do not expect a value.  Such keys still set 
                                                                         //some value themselves to be in the map of Policies.

      inline const std::string GetUsage() { return fUsage; }
      inline const std::string GetDetails() { return fDetails; }

      //Type of exception that will be thrown by all Policy-derived objects.  
      //TODO: Make a general exception template class that will take a type as a parameter.
      //      Then, each exception could be specific to the class that is throwing it without 
      //      rewriting exception code.       
      class exception: public util::exception
      {
        public:
          exception(const std::string& help) noexcept : util::exception("opt::Policy")
          {
            fExplanation << help;
          }
                                                                                          
          virtual ~exception() noexcept {};
      };

    protected:
      std::string fUsage; //Format of the type of argument that this option takes.
      std::string fDetails; //Long description of how to use this option.
  };
}

#endif //OPT_POL 
