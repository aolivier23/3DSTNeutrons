//File: GeoFunc.cpp
//Brief: Functions for dealing with ROOT's geometry facilities more conveniently. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoMatrix.h"

namespace geo
{
  TGeoMatrix* findMat(const std::string& name, TGeoNode& parent)
  {
    if(std::string(parent.GetVolume()->GetName()) == name) return parent.GetMatrix();
    auto children = parent.GetNodes();
    for(auto child: *children)
    {
      auto node = (TGeoNode*)child;
      auto result = findMat(name, *node);
      if(result)
      {
        auto retVal = new TGeoHMatrix(*(parent.GetMatrix()));
        retVal->Multiply(result);
        return retVal;
      }
    }
    return nullptr;
  }

  TVector3 InLocal(const TVector3& pos, TGeoMatrix* mat)
  {
    double master[3] = {}, local[3] = {};
    pos.GetXYZ(master);
    mat->MasterToLocal(master, local);
    return TVector3(local[0], local[1], local[2]);
  }

  TVector3 InGlobal(const TVector3& pos, TGeoMatrix* mat)
  {
    double master[3] = {}, local[3] = {};
    pos.GetXYZ(local);
    mat->LocalToMaster(local, master);
    return TVector3(master[0], master[1], master[2]);
  }

  double DistFromInside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    (begin-shapeCenter).GetXYZ(posArr);
    dir.GetXYZ(dirArr);

    return shape.DistFromInside(posArr, dirArr);
  }

  double DistFromOutside(const TGeoShape& shape, const TVector3& begin, const TVector3& end, const TVector3& shapeCenter)
  {
    const auto dir = (end-begin).Unit();
    double posArr[3] = {0}, dirArr[3] = {0};
    (begin-shapeCenter).GetXYZ(posArr);
    dir.GetXYZ(dirArr);

    return shape.DistFromOutside(posArr, dirArr);
  }

  bool Contains(const TGeoShape& shape, const TVector3& point, const TVector3& shapeCenter)
  {
    const auto diff = point-shapeCenter;
    double arr[] = {diff.X(), diff.Y(), diff.Z()};
    return shape.Contains(arr);
  }
} 
