//File: Options.cxx
//Brief: Access to command line options.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Options.h"

namespace opt
{
  Options::Options(const std::string& help, const std::string& exeName): fHelp(help), fExeName(exeName) {}

  /*std::string Options::operator [](const std::string& key) const
  {
    auto found = fOptions.find(key);
    if(found == cend()) throw exception(fHelp);
    return found->second;
  }

  template <class T>
  T Options::Get(const std::string& key) const
  {
    auto found = fOptions.find(key);
    if(found == cend()) throw exception(fHelp);
    return detail::decode<T>(found->second);
  }*/

  void Options::insert(const std::string& key, const std::string& value)
  {
    fOptions.insert(std::make_pair(key, value));
  }

  //TODO: This won't work with Counter policies.
  /*std::string Options::PrintToBash() const
  {
    std::string cmd = "#!/bin/bash\n"+fExeName;
    for(const auto& pair: fOptions)
    {
      cmd += " "+pair.first + " " + pair.second;
    }
    return cmd+"\n";
  }*/  
}
