//File: exception.h
//Brief: An exception interface for the util library.  Allows catching 
//       any exception from the util library by a reference to this 
//       base class.  All exception classes in util shall 
//       derive from this base class.  Inspired by the interface for 
//       cet::exception used in the ART framework.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//c++ includes
#include <sstream>

#ifndef UTIL_EXCEPTION_H
#define UTIL_EXCEPTION_H

//TODO: Finish making current exception classes derive from this one
namespace util
{
  class exception
  {
    public:
      exception(const std::string& category) noexcept;
      exception(const exception& other) noexcept;
      exception& operator =(const exception& rhs) noexcept;

      virtual ~exception() noexcept = default;

      virtual std::string what() const;
       
      template <class T>
      util::exception& operator << (const T& toPrint)
      {
        fExplanation << toPrint;
        return *this;
      }

    protected: 
      std::stringstream fExplanation; //Explaination of why this exception was thrown.
  };
}

#endif //UTIL_EXCEPTION_H
