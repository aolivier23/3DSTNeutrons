cmake_minimum_required( VERSION 3.5.1 )
project( EdepNeutrons )

#Compiler flags
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS ON )

set( GCC_Flags_For_CXX "-std=c++14 -Wall -fPIC" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )
set( CMAKE_CXX_FLAGS_RELEASE "-O2" )

#TODO: Does this cause problems if I relocate libraries like edepsim?  If I can just set LD_LIBRARY_PATH as a 
#      fallback, then put this back in.  Otherwise, remove RPATH stuff so that this package can be relocated to 
#      grid nodes.  
#Really!?!  Why is this a thing?
#Force CMake to link paths "correctly".  See the discussion at https://cmake.org/pipermail/cmake/2011-June/044807.html
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS Core Geom Tree RIO TreePlayer EG)
include(${ROOT_USE_FILE})

#Make sure that all project files know about all other project files by their 
#path from this directory.
include_directories( "${PROJECT_SOURCE_DIR}" )

add_subdirectory(persistency)
add_subdirectory(reco)
add_subdirectory(ana)
add_subdirectory(app)

#Make the results of this build into a package.  Designed to be distributed as a .tar.gz
#Learned to do this from http://agateau.com/2009/cmake-and-make-dist/
set( CPACK_PACKAGE_VERSION_MAJOR "1" )
set( CPACK_PACKAGE_VERSION_MINOR "0" )
set( CPACK_PACKAGE_VERSION_PATCH "0" )
set( CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
set( CPACK_SOURCE_GENERATOR "TGZ" )
#set( CPACK_SOURCE_IGNORE_FILES 

include( CPack )
add_custom_target( dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source )
