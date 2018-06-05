//File: TFileSentry.cxx
//Brief: A TFileSentry adds given TObjects to the TFile it manages.  It does not interfere with gFile. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "ROOT/Base/TFileSentry.h" 

namespace util
{
  TFileSentry::TFileSentry(const std::string& name): fFile(new TFile(name.c_str(), "RECREATE")) 
  {
  }

  TFileSentry::TFileSentry(std::shared_ptr<TFile> file): fFile(file) 
  {
  }

  TDirectory* TFileSentry::cd(const std::string& name)
  {
    auto oldDir = gDirectory;
    if(!fFile->cd(name.c_str()))
    {
      fFile->mkdir(name.c_str())->cd();
      fPwd = gDirectory; //TODO: There is probably a better way to extract a TDirectory from a TFile
      gDirectory = oldDir;
      return fPwd;
    }
    gDirectory = oldDir;
  }

  TFileSentry::~TFileSentry()
  {
    //Make sure all objects in this TFile are written.
    auto oldFile = gFile;
    fFile->cd(); 
    if(fFile->IsOpen() && fFile->GetList())
    {
      for(auto obj: *(fFile->GetList())) obj->Write();
    } 
    else throw util::exception("FileClosed") << "File " << fFile->GetName() << " was closed before TFileSentry was done "
                                             << "writing to it in destructor.\n";
    gFile = oldFile;
  }
}
