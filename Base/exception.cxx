//File: exception.cxx
//Brief: A base class for exceptions thrown by the util library.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "exception.h"

//c++ includes
#include <string>

namespace util
{
  exception::exception(const std::string& category) noexcept
  {
    fExplanation << category << ": ";
  }

  exception::exception(const exception& other) noexcept: fExplanation(other.fExplanation.str()) 
  {
  }

  exception& exception::operator =(const exception& rhs) noexcept
  {
    fExplanation.clear();
    fExplanation << rhs.what();
    return *this;
  }

  std::string exception::what() const
  {
    return fExplanation.str(); 
  }
}
