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
      #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
      if(traj.GetParentId() == parent)
      #else
      if(traj.ParentId == parent)
      #endif
      {
        #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
        const int id = traj.GetTrackId();
        #else
        const int id = traj.TrackId;
        #endif
        ids.insert(id);
        Descendants(id, trajs, ids);
      }
    }
  }
  
  //Returns the FS TG4Trajectory that led to child.
  const TG4Trajectory& Matriarch(const TG4Trajectory& child, const std::vector<TG4Trajectory>& trajs)
  {
    #ifdef EDEPSIM_FORCE_PRIVATE_FIELDS
    const int parent = child.GetParentId();
    #else
    const int parent = child.ParentId;
    #endif 
    if(parent == -1) return child;
    return Matriarch(trajs[parent], trajs);
  }
}
