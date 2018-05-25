//File: decode.cxx
//Brief: Standard string-to-anything conversion facility for util/IO library.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef DETAIL_DECODE_CXX
#define DETAIL_DECODE_CXX

//util includes
#include "Base/exception.h"

//c++ includes
#include <string>
#include <vector>
#include <algorithm>
 
namespace detail
{
  //Nasty template class that simplifies the code I have to write below.  Abandon all hope ye who enter here...
  template <typename TYPE>
  class Convertible
  {
    public:
      inline operator TYPE() const { return fValue; } 
    protected: 
      TYPE fValue;
  };

  //details for how to convert strings to user-requested type
  //Write your own overload for decode for types that do not fall into these categories
  //Exception thrown when decode() fails. 
  //TODO: Make this a member type when decode() is rewritten to become a class (so that I can support containers). 
  class decode_exception: public util::exception
  {
    public:
      decode_exception(const std::string& value, const std::string& msg) noexcept : exception("detail::decode")
      {
        fExplanation << "Could not decode value " << value << " because " << msg << "\n";
      }

      virtual ~decode_exception() noexcept {};

  };

  //Decode needs to be a class so I can use partial specializations to handle 
  //containers uniformly.  To not break the API I am already using for the 
  //non-class version, perform conversion in the constructor and provide a cast 
  //operator 
  template <typename TYPE>
  struct decode: public Convertible<TYPE>
  {
    decode(const std::string& value) { throw decode_exception(value, "No conversion exists for string "+value+" to requested type.\n"); }
  };

  template <>
  struct decode<std::string>: public Convertible<std::string>
  {
    decode(const std::string& value) { fValue = value; }
  };

  template <>
  struct decode<int>: public Convertible<int>
  {
    decode(const std::string& value)  
    {
      if(!std::all_of(value.begin(), value.end(), [](const char c) { return std::isdigit(c) || c == '+' || c == '-'; }))
      {
        throw decode_exception(value, "This string contains characters that are not numbers, +, or -.");
      }
      fValue = atoi(value.c_str());
    }
  };

  template<>
  struct decode<size_t>: public Convertible<size_t>
  {
    decode(const std::string& value)
    {
      if(!std::all_of(value.begin(), value.end(), [](const char c) { return std::isdigit(c) || c == '+' || c == '-'; }))
      {
        throw decode_exception(value, "This string contains characters that are not numbers, +, or -.");
      }
      fValue = atoi(value.c_str());
    }
  };

  template <>
  struct decode<double>: public Convertible<double>
  {
    decode(const std::string& value)
    {
      if(!std::all_of(value.begin(), value.end(), [](const char c) { return std::isdigit(c) || c == '+' 
                                                                         || c == '-' || c == 'e' || c == '.'; }))
      {
        throw decode_exception(value, "This string contains characters that are not numbers, +, -, ., or e (for scientific notation).");
      }
      fValue = atof(value.c_str());
    }
  };

  template <>
  struct decode<bool>: public Convertible<bool>
  {
    decode(const std::string& value)
    {
      if(value == "true" || value == "True" || value == "TRUE") fValue = true;
      else if(value == "false" || value == "False" || value == "FALSE") fValue = false;
      else 
      {
        try
        {
          fValue = decode<double>(value);
        }
        catch(const decode_exception& e)
        {
          throw decode_exception(value, "Got value "+value+" that could not be converted to bool.  Values that can be converted to "
                                       +"bool are: true, True, TRUE, false, False, FALSE, and any floating point number or integer.");
        }
      }
    }
  };

  //TODO: std::complex?  I've never used it, so I'm not yet convinced its needed.

  //Handle containers
  template<typename TYPE>
  struct decode<std::vector<TYPE>> //TODO: Why doesn't Convertible work here?  
  {
    public:
      decode(const std::string& value)
      {
        //This line defines the syntax for specifying a container on the command line, 
        //in a file, and anywhere else decode() is used.  
        const auto vals = tokenize(value);
        for(const auto& val: vals) fValue.push_back(decode<TYPE>(val));
      }

      inline operator std::vector<TYPE>() const
      {
        return fValue;
      }

    private:
      //TODO: Should I make this function available outside this class?  It seems like I rewrite it (incorrectly) pretty often...
      std::vector<std::string> tokenize(const std::string& value, const std::string& delim = " ")
      { 
        std::vector<std::string> tokens;
        size_t pos = 0, next = 0;
        while(next != std::string::npos)
        {  
          next = value.find(delim, pos);
          tokens.push_back(value.substr(pos, next-pos));
          pos = next+1;
        }
        return tokens;
      }

      std::vector<TYPE> fValue;
  }; 
}

#endif //DETAIL_DECODE_CXX
