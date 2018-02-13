//File: TruthFunc.cpp
//Brief: Shared functions for dealing with MC information.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include header
#include "alg/TruthFunc.h"

//edepsim includes
#include "TG4Trajectory.h"

namespace truth
{
  //TrackIDs of particles descended from parent returned by reference through ids.
  void Descendants(const int& parent, const std::vector<TG4Trajectory>& trajs, std::set<int>& ids)
  {
    for(const auto& traj: trajs)
    {
      if(traj.ParentId == parent)
      {
        const int id = traj.TrackId;
        ids.insert(id);
        Descendants(id, trajs, ids);
      }
    }
  }
  
  //Returns the FS TG4Trajectory that led to child.
  const TG4Trajectory& Matriarch(const TG4Trajectory& child, const std::vector<TG4Trajectory>& trajs)
  {
    if(child.ParentId == -1) return child;
    return Matriarch(trajs[child.ParentId], trajs);
  }
}
