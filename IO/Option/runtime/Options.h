//File: Options.h
//Brief: Options is mapping from key to value from the command line.  It is a wrapper over 
//       std::map<std::string, std::string> that throws on trying to access a key that is not 
//       present.  Exceptions thrown will print help text. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef OPT_OPTIONS_H
#define OPT_OPTIONS_H
 
//c++ includes
#include <map>
#include <string>

//util includes
#include "IO/Base/decode.cxx"
#include "Base/exception.h"

namespace opt
{
  class Options
  {
    public:
      Options(const std::string& help, const std::string& exeName);
      ~Options() = default;

      //direct element access
      std::string operator [](const std::string& key) const
      {
        auto found = fOptions.find(key);
        if(found == cend()) throw exception(fHelp, key);
        return found->second;
      }
  
      template <class T>
      T Get(const std::string& key) const
      {
        auto found = fOptions.find(key);
        if(found == cend()) throw exception(fHelp, key);
        return detail::decode<T>(found->second);
      }

      //iterator access
      /*using std::map<std::string, std::string>::iterator iterator;
      using std::map<std::string, std::string>::const_iterator const_iterator;*/

      std::map<std::string, std::string>::iterator begin() { return fOptions.begin(); }
      std::map<std::string, std::string>::iterator end() { return fOptions.end(); }
      std::map<std::string, std::string>::const_iterator cbegin() const { return fOptions.cbegin(); }
      std::map<std::string, std::string>::const_iterator cend() const { return fOptions.cend(); }

      std::map<std::string, std::string>::iterator find(const std::string& key) { return fOptions.find(key); }

      void insert(const std::string& key, const std::string& value);

      //std::string PrintToBash() const; //Prints a bash script to rerun this job with the same command line
 
      class exception: public util::exception
      {
        public:
          exception(const std::string& help, const std::string& key) noexcept : util::exception("opt::Options::exception")
          {
            fExplanation << "Key " << key << " was not specified.\n" << help;
          }
    
          virtual ~exception() noexcept {};
      };    

    private:
      std::map<std::string, std::string> fOptions;
      const std::string fHelp;
      const std::string fExeName;
  };
}

#endif //OPT_OPTIONS_H
