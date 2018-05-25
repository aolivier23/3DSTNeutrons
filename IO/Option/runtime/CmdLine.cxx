//File: CmdLine.cxx
//Brief: Implements a command line parser with user-definable Policies.  Access to options is provided through an external class, 
//       opt::Options, that is returned after parsing is complete.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//opt includes
#include "CmdLine.h"
#include "Policy.h"
#include "IO/Base/decode.cxx"
#include "Options.h"
#include "ExactlyOnce.h"
#include "Help.h"

//c++ includes
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace opt
{
  CmdLine::CmdLine(const std::string& desc): fDescription(desc) //Pass in a description of what this program does
  {
    //Add default help options.
    AddKey<Help>("--help", "");
    AddKey<Help>("-h", "");
  }

  CmdLine::~CmdLine() {}

  //TODO: Replace this with std::regex?
  bool CmdLine::isKey(const std::string& key) const
  {
    if(key.find('-') != 0) return false;
    return true;
    //TODO: What if key has 3 dashes?  Should I do anything about this?
  }

  opt::Options CmdLine::Parse(int argc, char** argv, const bool throwOnUnknown) //Parse data from the command line.   
  {
    const std::string exe(argv[0]);
    fUsage.insert(0, exe+" "); 
    opt::Options options(GetHelp(), exe); 
    int argCt = 1; //The first element of argv is always the name of the executable
    while(argCt < argc) //Loop over arguments
    {
      std::string key = argv[argCt];
      if(isKey(key))
      {
        auto foundPolicy = fPolicies.find(key);
        if(foundPolicy == fPolicies.end()) 
        {
          if(throwOnUnknown) throw exception("Got invalid key "+key+"\n"+GetHelp());
          else //Skip to the next key
          {
            ++argCt;
            if(argCt < argc && !isKey(argv[argCt])) ++argCt;
            continue; //TODO: Rewrite this to avoid using continue
          }
        }

        auto foundValue = options.find(key);

        //If the next token exists and is not itself a key, save it as the value for this key.
        ++argCt;
        std::string value("");
        if(argCt < argc && !isKey(argv[argCt])) 
        {
          value = std::string(argv[argCt]);
          ++argCt;

          try
          {
            if(foundValue == options.end()) options.insert(key, foundPolicy->second->FoundFirst(value));
            else 
            {
              //options[key] = foundPolicy->second->FoundAgain(value, foundValue->second);
              foundValue->second = foundPolicy->second->FoundAgain(value, foundValue->second);
            }
          }
          catch(const Policy::exception& e)
          {
            throw exception(std::string(e.what())+" for key "+key+"\n"+GetHelp()); 
          }
        }
        else //A key was specified without a value.  This is OK for some options.  Call the appropriate overload of FoundFirst/Again()
        {
          try
          {
            if(foundValue == options.end()) options.insert(key, foundPolicy->second->FoundFirst());
            else options[key] = foundPolicy->second->FoundAgain(foundValue->second);
          }
          catch(const Policy::exception& e)
          {
            throw exception(std::string(e.what())+" for key "+key+"\n"+GetHelp());
          }
        }  
        //Now, concatenate additional non-key values specified for this key.  
        //TODO: Perhaps I should require the user to put everything in '' like the rest of Unix?  
        //      This would substantially change the logic here.
        /*for(; argCt < argc && !isKey(argv[argCt]); ++argCt)
        {
          value += " "+std::string(argv[argCt]);
        }*/

        //Look for the policy to handle this key
        /*auto foundPolicy = fPolicies.find(key);
        if(foundPolicy == fPolicies.end()) throw exception("Got invalid key "+key+"\n"+GetHelp()); 

        auto foundValue = options.find(key);
        try
        {
          if(foundValue == options.end()) options.insert(key, foundPolicy->second->FoundFirst(value));
          else options[key] = foundPolicy->second->FoundAgain(value, foundValue->second);
        }
        catch(const Policy::exception& e)
        {
          throw exception(std::string(e.what())+" for key "+key+"\n"+GetHelp()); 
        }*/
      }
      else throw exception("Got invalid key "+key+"\n"+GetHelp()); 
    }

    ApplyDefaults(options);
    return options;
  }

  void CmdLine::ApplyDefaults(opt::Options& options) const
  {
    for(const auto& def: fDefaults)
    {
      if(options.find(def.first) == options.end()) options.insert(def.first, def.second);
    }
  }

  std::string CmdLine::GetHelp() const
  {
    std::stringstream stream;
    stream << "Usage:\n";
    stream << "   " << fUsage << "\n\n";
    stream << "Description:\n";
    stream << "   " << fDescription << "\n\n";
    stream << "Options: \n";    

    //Figure out the longest usage entry so I can insert appropriate spacing after the other entries.
    size_t maxLen = std::max_element(fHelp.begin(), fHelp.end(), 
                                     [](const std::pair<std::string, std::string>& lhs, const std::pair<std::string, std::string>& rhs) 
                                     { return lhs.first.size() < rhs.first.size(); })->first.size();

    //TODO: Move text wrapping to another function and wrap description text?
    //TODO: I do not think it would be difficult to add support for groups of options here.  Since the help text is the 
    //      only thing that needs to know about groups (other than the option-adding interface), I would just need to store 
    //      groups of options instead of option-description pairs here.  I should also probably add a simple option group class.  
    //      This would open up the opportunity to let Policy object reencode their arguments (Counter, I'm looking at you). 
    const size_t colLen = 60;
    for(auto& pair: fHelp)
    {
      stream << "   " << std::setw(maxLen+8) << std::left << pair.first;

      //Perform text wrapping manually.  There seems to be a Boost library to do this...
      auto desc = pair.second;

      if(colLen > desc.size()) stream << desc << "\n";
      else 
      {
        stream << desc.substr(0, desc.find(" ", colLen)) << "\n";
        size_t pos = desc.find(" ", colLen)+1;
        size_t nextPos = desc.find(" ", pos+colLen);
        while(nextPos < desc.size())
        {
          stream << std::setw(maxLen+11) << "" << desc.substr(pos, nextPos-pos) << "\n"; 
          pos = nextPos + 1; 
          nextPos = desc.find(" ", pos+colLen);
        }
        stream << std::setw(maxLen+11) << "" << desc.substr(pos) << "\n";
      }
    } 

    return stream.str();
  }
}
