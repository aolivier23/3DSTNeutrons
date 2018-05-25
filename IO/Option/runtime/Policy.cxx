//File: Policy.cxx
//Brief: A Policy handles how a command line option is handled.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//opt includes
#include "Policy.h"

namespace opt
{
  Policy::Policy(const std::string& usage, const std::string& details): fUsage(usage), fDetails(details) {}
  //Trivial constructor to simultaneously make CMake happy and maintain my file include convention.
}
