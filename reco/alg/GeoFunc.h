//File: GeoFunc.h
//Brief: Shared functions for answering the questions I am 
//       actually inerested in with ROOT's TGeoManager and friends.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef GEO_GEOFUNC_H
#define GEO_GEOFUNC_H

//Forward declarations for ROOT classes
class TGeoNode;
class TVector3;
class TLorentzVector;
class TGeoMatrix;
class TGeoShape;

//c++ includes
#include <string>

namespace geo
{
  //Return the product of the matrices from the node with a volume called name with all of its' ancestors.
  TGeoMatrix* findMat(const std::string& name, TGeoNode& parent);
                                                                                                                          
  //Convert a 3-vector from a global position to the local coordinate system described by mat
  TVector3 InLocal(const TVector3& pos, TGeoMatrix* mat);
                                                                                                                          
  //The opposite of InLocal
  TVector3 InGlobal(const TVector3& pos, TGeoMatrix* mat);
                                                                
  //Finds the distance until leaving a boundary of shape for a line from begin to end.  Begin must be inside shape.                 
  double DistFromInside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter);
                                                                                                                  
  //Finds the distance until entering a boundary of shape for a line from begin to end.  Begin must be outside shape.        
  double DistFromOutside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter);
  
  //Checks whether a point is inside shape.
  bool Contains(const TGeoShape& shape, const TVector3& point, const TVector3& shapeCenter);
}
#endif //GEO_GEOFUNC_H
