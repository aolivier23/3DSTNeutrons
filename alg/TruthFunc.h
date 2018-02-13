//File: TruthFunc.h
//Brief: Shared functions for dealing with truth information.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//c++ includes
#include <vector>
#include <set>

class TG4Trajectory;

namespace truth
{
  //Get a set of the particles that are descended from parent.
  void Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::set<int>& ids);

  //Get the FS trajectory that created child.  
  const TG4Trajectory& Matriarch(const TG4Trajectory& child, const std::vector<TG4Trajectory>& trajs);
}
