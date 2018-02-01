//File: SpaceSort.cpp
//Brief: Algorithm to break up objects with 3-vectors into octants.  Makes between 3*(subdiv +1)
//       decisions per object to sort objects into 2^(subdiv*3) octants.  Access an object's assignment in 
//       3*(subdiv +1) decisions.  So, search is linear in number of queries and constant in 
//       number of objects already sorted.  Roughly linear in number of subdivisions.  With lines like TG4HitSegment, 
//       just use their centers.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TVector3.h"

namespace
{
  template <class VEC, size_t COMP>
  TVector3 Add(VEC& vec, const double width)
  {
    vec[COMP] += width;
    return vec;
  }
}

namespace reco
{
  //User interface
  template <class CELL, size_t DEPTH, class KEY=TVector3>
  class OcTree
  {
    public:
      OcTree(const KEY& center, const double width): fCenter(center), fWidth(width) {}

      //CELL* is an observer pointer.  
      std::pair<KEY, CELL*> operator [](const KEY& pos)
      {
        auto pos = fCenter;
        auto cell = fRoot[pos, pos, fWidth]; 
        return std::make_pair(pos, cell);
      }

      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        fRoot.visitor(visit);
      }
    
    private:
      KEY fCenter;
      double fWidth;
      Node<CELL, 0, DEPTH> fRoot;
  }; 

  //Compile-time implementation.  Very compact, but depth is fixed at compile-time.
  template <class CELL, size_t COMP, size_t LEVEL, class KEY=TVector3>
  class Node
  {
    private:
      Node<CELL, COMP+1, LEVEL> Plus;
      Node<CELL, COMP+1, LEVEL> Minus;
    
    public:  
      CELL* operator[](const KEY& pos, KEY& center, const double width)
      {
        return (pos[COMP] < center[COMP])?Minus[pos, ::Add<COMP>(center, -width/2.)]:Plus[pos, ::Add<COMP>(center, width/2.)];
      }
    
      //TODO: Can I use SFINAE to allow this to return a value based on VISITOR's return type?
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        Plus.visitor(visit);
        Minus.visitor(visit);
      }
  };
  
  //Partial specialization for non-terminal z node
  template <class CELL, size_t, size_t LEVEL, class KEY=TVector3>
  class Node<CELL, 2, LEVEL>
  {
    private:
      Node<CELL, 0, LEVEL-1> Plus;
      Node<CELL, 0, LEVEL-1> Minus;
    
    public:  
      CELL* operator [](const KEY& pos, KEY& center, const double width)
      {
        return (pos[2] < Center[2])?Minus[pos, ::Add<2>(center, -width/2.)]:Plus[pos, ::Add<2>(center, width/2.)];
      }
    
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        Plus.visitor(visit);
        Minus.visitor(visit);
      }
  };
  
  //Partial specialization for terminal z node
  template <class CELL, size_t, size_t, class KEY=TVector3>
  class Node<CELL, 2, 0>
  {
    private:
      std::unique_ptr<CELL> Plus;
      std::unique_ptr<CELL> Minus;
    
      KEY Center;
  
    public: 
      CELL* operator[](const KEY& pos, KEY& /*center*/, const double /*width*/)
      {
        return (pos[2] < Center[2])?Minus.get():Plus.get();
      }
  
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        visit(Plus);
        visit(Minus);
      }
  };
}
