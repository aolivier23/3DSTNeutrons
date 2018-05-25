//File: TCollectionSTLIter.h
//Brief: A wrapper over TIter to provide an STL-style iterator over a TCollection.  Should also work with any TCollection-derived class 
//       (i.e. all of the containers in ROOT).
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOT_TCOLLECTIONSTLITER_H
#define PLOT_TCOLLECTIONSTLITER_H

//c++ includes
#include <memory>

//Forward declarations of some ROOT classes
class TObject;
class TCollection;
class TIter;

namespace plot
{
  class TCollectionSTLIter
  { 
    public:
      TCollectionSTLIter(TCollection& list);
      TCollectionSTLIter(TIter& iter);
      TCollectionSTLIter(TIter&& iter);
      TCollectionSTLIter(TCollectionSTLIter& other); //copy constructor

      ~TCollectionSTLIter();

      //Copy facilities
      TCollectionSTLIter& operator =(TCollectionSTLIter& other);
      
      //Comparison
      inline bool operator ==(const TCollectionSTLIter& other) const { return other.fNext == fNext; }
      inline bool operator !=(const TCollectionSTLIter& other) const { return !(other == *this); }

      //Increment and decrement(?)
      TCollectionSTLIter operator ++();
      TCollectionSTLIter operator ++(int);

      //Object access
      TObject* operator *(); //Returning TObject* because I need to be able to cast back to the object I want to work with (unfortunately).  
                             //This make TCollectionSTLIter "look" like an iterator to a container of TObject pointers.  

    private:
      std::unique_ptr<TIter> fNext;
  };

   
  class const_TCollectionSTLIter
  { 
    public:
      const_TCollectionSTLIter(const TCollection& list);
      const_TCollectionSTLIter(const TIter& iter);
      const_TCollectionSTLIter(const TIter&& iter);
      const_TCollectionSTLIter(const const_TCollectionSTLIter& other); //copy constructor
                                                                                                                                              
      ~const_TCollectionSTLIter();
                                                                                                                                              
      //Copy facilities
      const_TCollectionSTLIter& operator =(const const_TCollectionSTLIter& other);
      
      //Comparison
      inline bool operator ==(const const_TCollectionSTLIter& other) const { return other.fNext == fNext; }
      inline bool operator !=(const const_TCollectionSTLIter& other) const { return !(other == *this); }
                                                                                                                                              
      //Increment and decrement(?)
      const_TCollectionSTLIter operator ++();
      const_TCollectionSTLIter operator ++(int);
                                                                                                                                              
      //Object access
      const TObject* operator *(); //Returning TObject* because I need to be able to cast back to the object I want to work with (unfortunately).  
                                   //This make TCollectionSTLIter "look" like an iterator to a container of TObject pointers.  
                                                                                                                                              
    private:
      std::unique_ptr<TIter> fNext; //Not using unique_ptr because I want to be able to copy this class
  };

}


//Since TCollection doesn't work well with the STL, the following free functions are needed to support range-based for loops over TCollection using TCollectionSTLIter
plot::TCollectionSTLIter begin(TCollection& list);
plot::TCollectionSTLIter end(TCollection& list);

plot::const_TCollectionSTLIter cbegin(const TCollection& list);
plot::const_TCollectionSTLIter cend(const TCollection& list);

#endif //PLOT_TCOLLECTIONSTLITER_H
