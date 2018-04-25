//File: NeutronCand.h
//Brief: A NeutronCand is a group of MCClusters that I hypothesize are from the same (FS?) neutron.  It has a guess at neutron energy from Time Of Flight, the sum of 
//       deposited energy from its' clusters, the indices of the clusters it contains for each clustering algorithm used, and the list of unique true particles that 
//       contributed to its' clusters.  I might add a direction eventually, especially if this is a candidate for a FS neutron.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include <TObject.h>
#include <TLorentzVector.h>

//c++ includes
#include <set>
#include <map>

#ifndef PERS_NEUTRONCAND_H
#define PERS_NEUTRONCAND_H

namespace pers
{
  class NeutronCand: public TObject
  { 
    public:
      NeutronCand(void) {}
      virtual ~NeutronCand(); //This needs to be defined in the .cpp file so that I can call ClassImp(?)
      
      //Energy deposited in this NeutronCand
      double TOFEnergy; //Energy estimated from TOF of clusters
      double Beta; //Velocity from TOF and distance from vertex over speed of light 
      double SigmaBeta; //Uncertiainty on Beta
      double DepositedEnergy; //Energy deposited by clusters
      
      //"Beginning" of this NeutronCand as defined by algorithm that produced it.  Will be used to study angle of candidate w.r.t. true neutron.
      TLorentzVector Start; //All in mm like edep-sim

      //TrackIDs that contributed to this NeutronCand.  
      std::set<int> TrackIDs; //TODO: Do I need this if I can backtrack to MCClusters?   

      //Indices of MCClusters in this neutron candidate
      std::map<std::string, std::vector<size_t>> ClusterAlgToIndices; //TODO: Turn these indices into clusters/links/iterators?  I could do this like what I think ART 
                                                                      //      does and store everything in a Ptr class template that has a unique index, but I'd rather 
                                                                      //      my file format stay as close to edep-sim as possible.  

      ClassDef(NeutronCand, 1);
  };
}

#endif //PERS_NEUTRONCAND_H
