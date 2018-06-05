//File: TFileSentry.h
//Brief: A TFileSentry makes sure I know what TFile objects I want to write to a ROOT file are actually 
//       written to.  It should be careful not to interfere with other TFiles. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TFile.h"

//util includes
#include "Base/exception.h"

//c++ includes
#include <memory>
 
#ifndef UTIL_TFILESENTRY_H
#define UTIL_TFILESENTRY_H

//This concept is based on LArSoft's TFileService.  
namespace util
{
  class TFileSentry
  {
    public:
      //For convenience, create a new TFile shared_ptr.  fFile is still publicly accessible.
      TFileSentry(const std::string& name);
                                                                                             
      //Create from an existing std::shared_ptr<TFile>.
      TFileSentry(std::shared_ptr<TFile> file);
                                                                                             
      template <class TOBJECT, class ...ARGS>
      TOBJECT* make(ARGS... args)
      {
        auto oldFile = gFile;
        auto oldDir = gDirectory;
        fPwd->cd();
        if(fFile->IsOpen())
        {
          auto obj = new TOBJECT(args...);
          gFile = oldFile;
          gDirectory = oldDir;
          return obj;
        }
        else throw util::exception("FileClosed") << "File " << fFile->GetName() << " was closed before TFileSentry was done "
                                                 << "writing to it in make.\n";
        gFile = oldFile;
        gDirectory = oldDir;
        return nullptr;
      }

      //Change the directory for newly written TObjects the the TDirectory called "name", creating it if necessary.
      TDirectory* cd(const std::string& name); 
                                                                                             
      virtual ~TFileSentry();
  
      std::shared_ptr<TFile> fFile; //The file to which objects will be written.  
      TDirectory* fPwd; //Observer pointer to current working directory.  
  };
}

#endif //UTIL_TFILESENTRY_H
