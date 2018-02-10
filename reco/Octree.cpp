//File: SpaceSort.cpp
//Brief: Algorithm to break up objects with 3-vectors into octants.  Makes between 3*(subdiv +1)
//       decisions per object to sort objects into 2^(subdiv*3) octants.  Access an object's assignment in 
//       3*(subdiv +1) decisions.  So, search is linear in number of queries and constant in 
//       number of objects already sorted.  Roughly linear in number of subdivisions.  With lines like TG4HitSegment, 
//       just use their centers.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TVector3.h"

//c++ includes
#include <memory>

namespace reco
{
  //Compile-time implementation.  Very compact, but depth is fixed at compile-time.
  template <class CELL, size_t COMP, size_t LEVEL, class KEY=TVector3>
  class Node
  {
    private:
      Node<CELL, COMP+1, LEVEL> Plus;
      Node<CELL, COMP+1, LEVEL> Minus;
    
    public:  
      CELL* Get(const KEY& pos, KEY& center, const KEY width)
      {
        return (pos[COMP] < center[COMP])?Minus.Get(pos, center, width):Plus.Get(pos, center, width);
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
  template <class CELL, size_t LEVEL, class KEY>
  class Node<CELL, 2, LEVEL, KEY>
  {
    private:
      Node<CELL, 0, LEVEL-1> Plus;
      Node<CELL, 0, LEVEL-1> Minus;
    
    public:  
      CELL* Get(const KEY& pos, KEY& center, const KEY width)
      {
        if(pos[2] < center[2])
        {
          center += width*(-0.5);
          return Minus.Get(pos, center, width*0.5);
        }

        //else
        center += width*0.5;
        return Plus.Get(pos, center, width*0.5); //Update width for next level of Nodes
      }
    
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        Plus.visitor(visit);
        Minus.visitor(visit);
      }
  };
  
  //Partial specialization for terminal z node
  template <class CELL, class KEY>
  class Node<CELL, 2, 0, KEY>
  {
    private:
      std::unique_ptr<CELL> Plus;
      std::unique_ptr<CELL> Minus;
    
    public: 
      Node(): Plus(nullptr), Minus(nullptr) {}

      CELL* Get(const KEY& pos, KEY& center, const KEY /*width*/)
      {
        return (pos[2] < center[2])?Minus.get():Plus.get();
      }
  
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        visit(Plus.get());
        visit(Minus.get());
      }
  };

  //User interface
  template <class CELL, size_t DEPTH, class KEY=TVector3>
  class Octree
  {
    public:
      Octree(const KEY& center, const KEY& width): fCenter(center), fWidth(width), fRoot() {}
                                                                                                                       
      //CELL* is an observer pointer.  
      std::pair<KEY, CELL*> operator [](const KEY& pos) //Make this look somewhat like a std::map even though it isn't
      {
        auto center = fCenter;
        auto cell = fRoot.Get(pos, center, fWidth); 
        return std::make_pair(center, cell);
      }
                                                                                                                       
      //This is how you should loop over the elements of an Octree (unfortunately).  
      //VISITOR is any callable object, so it can be stateful.  
      //You could also use a lambda function that just captures the variable for 
      //the result you want to calculate.  
      template <class VISITOR>
      void visitor(VISITOR& visit)
      {
        fRoot.visitor(visit);
      }

      template <class VISITOR>
      void visitor(VISITOR&& visit)
      {
        fRoot.visitor(visit);
      }
    
    private:
      const KEY fCenter;
      const KEY fWidth;
      Node<CELL, 0, DEPTH> fRoot;
  }; 
}
