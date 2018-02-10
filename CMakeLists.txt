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

#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS Core Geom Tree RIO TreePlayer EG)
include(${ROOT_USE_FILE})

#Variable for location of the util package
#TODO: Make util a real package so I don't need to do this
#set( UTIL_DIR "/dune/app/users/aolivier/nd_neutrons/util")
message("Linking util library at ${UTIL_DIR}")

#Make sure output goes to the right place
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../lib )
set( CMAKE_INSTALL_LIBDIR           ${PROJECT_SOURCE_DIR}/../lib )

#Make sure that all project files know about all other project files by their 
#path from this directory.
include_directories( "${PROJECT_SOURCE_DIR}" )

add_subdirectory(persistency)
add_subdirectory(reco)
add_subdirectory(ana)
add_subdirectory(app)
