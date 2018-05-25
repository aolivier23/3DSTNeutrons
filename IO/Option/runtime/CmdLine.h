//File: CmdLine.h
//Brief: A CmdLine parses the command line and stores the options it finds in its own mapping.  
//       
//       Command line behavior:
//         In general, CmdLine expects options to be specified with strings like
//         --option, --option-name or, -option.  The string after --option-name and before 
//         the next --another_option will be given to the option object for --option-name 
//         to be parsed.  If an option gets input that it cannot parse, an exception will 
//         be thrown by CmdLine that will end the program with help information unless 
//         user code catches it.    
//
//       Element Access:
//         There are three ways to ask for options from the command line.  To access options directly, 
//         use the access operator, [], or Get().  Both of these functions will throw a CmdLine::exception 
//         to print usage information if the requested key was not found.  There are also begin() and 
//         end() methods which return iterators to the map of keys to options.  
//
//       Configuration:
//         Keys to look for are added using the method AddKey().  AddKey() supports user-specified 
//         option types via the Policy class provided with this package.  There are several AddKey() 
//         overloads that support a few basic option types.  The default AddKey(const std::string& key)
//         adds a key associated with an option that must be specified exactly once. 
// 
//         template <class POL=opt::ExactlyOnce, class ...ARGS>
//         AddKey<OPT>(const std::string& key, const std::string& default, ARGS... args)
//
//         sets a default value for the provided key AFTER parsing in complete if the key was never found. args are 
//         used to construct a Policy-derived object.
//
//       Policy:
//         Policy-derived objects use two methods to determine how to handle command line options for the keys 
//         they are paired with.  FoundFirst() is called the first time a key is found.  FoundAgain() is called 
//         the second time a key is found.  I have provided a few helpful Policy-derived classes in this package:
//         
//         Policy:      The base class for Policy classes.  It is an abstract class, so it is only an interface.  
//                      A Policy contains help information.
//         ExactlyOnce: Expects a string to be provided from the command line exactly once.  Will throw in FoundAgain().  
//                      Throwing if no default was set at the end of parsing is handled by operator[]() and Get().
//         Counter:     Expects an option only with an empty string provided.  Counts the number of times its key 
//                      appears.
//         Accumulate:  Expects a string to be provided from the command line each time its key appears.  Concatenates 
//                      strings with a space between them in FoundAgain().
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_CMDLINE_H
#define OPT_CMDLINE_H

//util includes
#include "Base/exception.h"

//c++ includes
#include <string>
#include <set>
#include <memory>
#include <map>
#include <vector>

namespace opt
{
  class Policy;
  class Options; 
  class ExactlyOnce;

  class CmdLine
  {
    public:
      CmdLine(const std::string& desc);
      virtual ~CmdLine(); //c++ quirk: since I am using smart pointers, I need to define the destructor in a scope where 
                          //this class's source code is already linked to all objects it uses.  That usually means the 
                          //.cxx file.   

      //The type of exception thrown by this object.  Any opt::Policy::exceptions caught 
      //will be printed before throwing this exception.
      class exception: public util::exception
      {
        public:
          exception(const std::string& help) noexcept : util::exception("opt::CmdLine") 
          {
            fExplanation << help;
          }

          virtual ~exception() noexcept {};
      };
      
      opt::Options Parse(int argc, char** argv, const bool throwOnUnknown = true); //Parse data from the command line.   Return a lightweight 
                                                                                   //object with the final command line options.  You should 
                                                                                   //probably destroy this object after obtaining the opt::Options 
                                                                                   //from Parse.

      //Configuration
      template <class POL=ExactlyOnce>
      void AddKey(const std::string& key, const std::string& desc) //Add an option with no default value.  desc explains what this key's 
                                                                   //option is used for.  Usage syntax is included included for you.
      {
        auto policy = new POL;
        fPolicies.insert(std::make_pair(key, std::unique_ptr<Policy>(policy)));
        const std::string usage = key+" "+policy->GetUsage();
        fUsage += usage+" ";
        //fDetails += usage+policy->GetDetails()+desc+"\n"; //TODO: Add correct number of spaces to help message. 
        fHelp.push_back(make_pair(usage, policy->GetDetails()+desc));
      }

      template <class POL=ExactlyOnce>
      void AddKey(const std::string& key, const std::string& desc, const std::string& def) //Set a default value for this option
      {
        auto policy = new POL;
        fPolicies.insert(std::make_pair(key, std::unique_ptr<Policy>(policy)));
        const std::string usage = key+" "+policy->GetUsage()+" "+"[="+def+"]";
        fUsage += usage+" "; //Include default value in usage information
        //fDetails += usage+": "+policy->GetDetails()+desc+"\n"; //TODO: Add correct number of spaces to help message.
        fHelp.push_back(std::make_pair(usage, policy->GetDetails()+desc));
        fDefaults.insert(std::make_pair(key, def));
      }

      std::string GetHelp() const;

      //Iterators
      /*using std::map<std::string, std::string>::iterator;
      using std::map<std::string, std::string>::const_iterator;

      iterator begin() const { return fOptions.begin(); }
      iterator end() const { return fOptions.end(); }

      const_iterator cbegin() const { return fOptions.cbegin(); }
      const_iterator cend() const { return fOptions.cend(); }

      //Direct element access
      std::string& operator [](const std::string& key) const;
      
      template <class T>
      T& Get(const std::string& key); */

    private:
      void ApplyDefaults(Options& options) const; //Apply default values just before returning from Parse. 
      bool isKey(const std::string& arg) const; //Check if arg is a valid key.  Enforce syntax rules of --key, --key-name, and -key here.

      std::map<std::string, std::string> fDefaults; //Default values
      std::map<std::string, std::unique_ptr<Policy>> fPolicies; //Map of keys to Policies      
      std::string fUsage; //Lists the options that can be specified for this program
      std::string fDescription; //A description of what this option parser does
      //std::string fDetails; //Gives details for each option
      std::vector<std::pair<std::string, std::string>> fHelp; //Pair of (usage, description) pairs.
      //std::map<std::string, std::string> fOptions; //Map of key to option from the command line.
      //TODO: Implement option to print command line to a file. 
  };
}

#endif //OPT_CMDLINE_H
