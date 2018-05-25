//File: Help.cxx
//Brief: A Help option throws an exception any time its key is found.  This will trigger printing help information.  
//       Intended for automatic generation of -h and --help options, but could be used by the user for other things 
//       as well.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "Help.h"

namespace opt
{
  Help::Help(): Policy("", "Prints this help text.")
  {
  }

  std::string Help::FoundFirst(const std::string& cmdLine) const 
  {
    throw exception("Help requested"); //Trigger printing of help text by CmdLine class
  }

  std::string Help::FoundAgain(const std::string& cmdLine, const std::string& prev) const
  {
    throw exception("Help requested"); //Trigger printing of help text by CmdLine class
  }

  std::string Help::FoundFirst() const
  {
    throw exception("Help requested"); //Trigger printing of help text by CmdLine class
  }

  std::string Help::FoundAgain(const std::string& prev) const
  {
    throw exception("Help requested"); //Trigger printing of help text by CmdLine class
  }
}
