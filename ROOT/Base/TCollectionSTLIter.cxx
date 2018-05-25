//File: TCollectionIter.cpp
//Brief: A wrapper over TIter to make an STL-style iterator over a TCollection.  "Looks" like an iterator over a container of TObjects.  Combined with 
//       a few overloads of std namespace functions, this allows for range-based for loops over TCollections.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

//local includes
#include "TCollectionSTLIter.h"

//ROOT includes
#include "TCollection.h"
#include "TObject.h"
#include "TCollection.h"

namespace plot
{
  TCollectionSTLIter::TCollectionSTLIter(TCollection& list): fNext(new TIter(&list))
  {
  }

  TCollectionSTLIter::TCollectionSTLIter(TIter& iter): fNext( new TIter(iter))
  {
  }

  TCollectionSTLIter::TCollectionSTLIter(TIter&& iter): fNext( new TIter(iter))
  {
  }

  TCollectionSTLIter::TCollectionSTLIter(TCollectionSTLIter& other): fNext(other.fNext.get())
  {
  }
  
  TCollectionSTLIter::~TCollectionSTLIter() {}

  TCollectionSTLIter& TCollectionSTLIter::operator =(TCollectionSTLIter& other)
  {
    fNext.reset(other.fNext.get());
    return *this;
  }

  TCollectionSTLIter TCollectionSTLIter::operator ++()
  {
    (*fNext)();
    return *this;
  }

  TCollectionSTLIter TCollectionSTLIter::operator ++(int)
  {
    (*fNext)();
    return *this;
  }

  TObject* TCollectionSTLIter::operator *()
  {
    return *(*fNext);
  }

  
  const_TCollectionSTLIter::const_TCollectionSTLIter(const TCollection& list): fNext(new TIter(&list))
  {
  }

  const_TCollectionSTLIter::const_TCollectionSTLIter(const TIter& iter): fNext(new TIter(iter))
  {
  }
  
  const_TCollectionSTLIter::const_TCollectionSTLIter(const TIter&& iter): fNext(new TIter(iter))
  {
  }

  const_TCollectionSTLIter::const_TCollectionSTLIter(const const_TCollectionSTLIter& other): fNext(other.fNext.get())
  {
  }

  const_TCollectionSTLIter::~const_TCollectionSTLIter() {}

  const_TCollectionSTLIter& const_TCollectionSTLIter::operator =(const const_TCollectionSTLIter& other)
  {
    fNext.reset(other.fNext.get());
    return *this;
  }

  const_TCollectionSTLIter const_TCollectionSTLIter::operator ++()
  {
    (*fNext)();
    return *this;
  }

  const_TCollectionSTLIter const_TCollectionSTLIter::operator ++(int)
  {
    (*fNext)();
    return *this;
  }

  const TObject*   const_TCollectionSTLIter::operator *()
  {
    return *(*fNext);
  }
}


//Since TCollection doesn't work well with the STL, the following free functions are needed to support range-based for loops over TCollection using TCollectionSTLIter
plot::TCollectionSTLIter begin(TCollection& list)
{
  plot::TCollectionSTLIter it(list);
  return it;
}

plot::TCollectionSTLIter end(TCollection& /*list*/)
{
  plot::TCollectionSTLIter it(TIter::End());
  return it;
}

plot::const_TCollectionSTLIter cbegin(const TCollection& list)
{
  plot::const_TCollectionSTLIter it(list);
  return it;
}

plot::const_TCollectionSTLIter cend(const TCollection& /*list*/)
{
  plot::const_TCollectionSTLIter it(TIter::End());
  return it;
}
