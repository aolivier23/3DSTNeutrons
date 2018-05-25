//File: RegexFiles.cxx
//Brief: Make a collection of arbitrary file objects based on a regular expression.  
//Author: Andrew Olivier aolivier@ur.rochester.edu
//Date: 1/24/2017

//util includes
#include "ROOT/Base/TCollectionSTLIter.h"

//c++ includes
#include <regex>
#include <string> //for std::to_string?
#include <vector>
#include <iostream> 

//ROOT includes
#include "TSystemDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TSystem.h"

#ifndef REGEXFILES_H
#define REGEXFILES_H

//Nasty details coming in the anonymous namespace!  
//Basically, this makes sure that vectors of pointer 
//types use new in emplace_back.
namespace
{
  //A FileAppender decides whether to create 
  //a file or call new to create a pointer 
  //to a file based on whether T is a pointer.
  template <class T, class ...ARGS>
  struct FileAppender
  {
    static void Append(std::vector<T>& vec, const std::string& name, ARGS... args)
    {
      vec.emplace_back(name, args...);
    }
  };

  //partial specialization for pointer types
  template <class T, class ...ARGS>
  struct FileAppender<T*, ARGS...>
  {
    static void Append(std::vector<T>& vec, const std::string& name, ARGS... args) //TODO: This shouldn't work!  My guess is that it has never been tested
    {
      vec.emplace_back(new T(name, args...));
    }
  };
}

namespace util
{
  template <class FILETYPE=std::string, class ...ARGS> 
  std::vector<FILETYPE> RegexFiles(TSystemDirectory& pwd, const std::regex& regex, std::vector<FILETYPE>& files, 
                                   const bool matchPath = false, ARGS... args)
  {
    std::string path = std::string(pwd.GetTitle())+"/";
    //Each TSystemFile knows its path (its title), so finding path here isn't strictly needed.  
    //However, I would have to add an additional function call to get each file's 
    //path, and this function would return exactly what path is currently set to.  
    //So, calculate path only once here to reduce the work to be done in the loop 
    //over files.
    //std::cout << "Path is " << path << "\n";
    for(auto obj: *(pwd.GetListOfFiles())) 
    {
      auto file = (TSystemFile*)(obj); 
      if(std::string(file->GetName()) == "." || std::string(file->GetName()) == "..") continue; //Ignore the current and parent directories
      //TODO: Second-guess using continue
      std::string fName = path+std::string(file->GetName()); 
      //std::cout << "In RegexFiles, examining file named " << fName << "\n";
      //TODO: It would be useful to have some kind of partial regex here before descending to a directory.
      if(file->IsDirectory()) RegexFiles(*((TSystemDirectory*)(file)), regex, files, matchPath, args...); //Search this directory recursively
      else //This is a file, not a directory.  So, determine whether it matches regex.
      {
        if(matchPath && (std::regex_match(fName, regex))) FileAppender<FILETYPE>::Append(files, fName, args...); //Match full path
        else if(std::regex_match(file->GetName(), regex)) FileAppender<FILETYPE>::Append(files, fName, args...); //Match file name only
        //else std::cout << "File " << fName << " failed to match regular expression \n"; 
      }
    }
    //std::cout << "In RegexFiles, the size of the vector of files to be returned is " << files.size() << "\n";
    return files; //TODO: I am returning files both by value and by reference!  Decide which one to use.  By reference 
                  //      seems more memory efficient for recursion. 
  }

  //TODO: It might be smarter to have one overload that checks for absolute/relative paths and handles them appropriately.
  template <class FILETYPE=std::string, class ...ARGS> 
  std::vector<FILETYPE> RegexFilesPath(const std::string& regex, const std::string& path, const bool matchPath = false, ARGS... args)
  {
    //std::cout << "Comparing to regular expression " << regex << "\n";
    std::vector<FILETYPE> files;

    //Separate the directory I am interested in and its path
    const size_t lastSlash = path.find_last_of('/');
    const auto dirName = path.substr(lastSlash+1, std::string::npos); //Directory name by itself

    TSystemDirectory dir(dirName.c_str(), path.c_str()); 
    /*std::cout << "The files in RegexFilesPath's starting directory, " << dir.GetTitle() << " are:\n";
    auto fileList = dir.GetListOfFiles();
    if(fileList == nullptr) std::cout << "Nothing!  There was an error in TSystemDirectory that returned a nullptr!.\n";
    else fileList->Print();*/
    //std::cout << "In RegexFiles, about to start matching to regular expression " << regex << "\n";
    RegexFiles(dir, std::regex(regex), files, matchPath, args...); 
  
    //std::cout << "In RegexFilesPath, the size of the vector of files to be returned is " << files.size() << "\n";
    return files;
  }

  template <class FILETYPE=std::string, class ...ARGS>
  std::vector<FILETYPE> RegexFiles(const std::string& regex, const bool matchPath = false, ARGS... args)
  {
    std::vector<FILETYPE> files;
    const auto path = std::string(gSystem->pwd()); 
    return RegexFilesPath((matchPath?path+"/":"")+regex, path, matchPath, args...); 
  }
}

#endif //REGEXFILES_H
