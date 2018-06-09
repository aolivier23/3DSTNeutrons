cmake_minimum_required( VERSION 3.5.1 )
project( 3DSTNeutrons )

#Compiler flags
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS ON )

set( GCC_Flags_For_CXX "-std=c++14 -Wall -fPIC" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )
set( CMAKE_CXX_FLAGS_RELEASE "-O2" )

#Find the EDepSim library that this package is designed to work with.  I will provide for 3 different strategies:
#1.) Look for the currently setup edep-sim.  This is not CMake's usual/default behavior as far as I can tell.  
#2.) Look for edep-sim that is "registered" via CMake.  I think this is how CMake usually looks for packages.  
#3.) Look for edep-sim without a CMake setup system.  v2.0.1 and older didn't have a Find*.cmake file.  

#First, look for currently setup edep-sim.  This should override all other cases/
if(DEFINED ENV{"$ENV{EDEP_ROOT}/$ENV{EDEP_TARGET}/lib"})
  find_package(EDepSim HINTS "$ENV{EDEP_ROOT}/$ENV{EDEP_TARGET}/lib" NO_DEFAULT_PATH)
else()
  find_package(EDepSim) #If edep-sim is not set up, try to find it "the CMake way" as documented at 
                        #https://cmake.org/cmake/help/v3.0/command/find_package.html
endif()

if(EDepSim_FOUND) #If edep-sim is not new enough to have CMake support, set it up using EDEP_DIR
  add_definitions(-DEDEPSIM_FORCE_PRIVATE_FIELDS=true) #Use accessor functions instead of public fields so that Clark can change the format of 
                                                       #the io classes without breaking things like edepViewer
  set(EDepSimIO EDepSim::edepsim_io)
  include_directories("${EDepSim_DIR}/include/EDepSim")
  message("EDepSim_DIR is " ${EDepSim_DIR} )
else() #If I couldn't find edep-sim using either its' environment variables or CMake's default behavior, try the old edep-sim setup method
  #add_definitions(-DEDEPSIM_FORCE_PRIVATE_FIELDS=false) #This option wasn't available until after v2.0.1
  set(EDepSimIO edepsim)
  link_directories("$ENV{EDEP_ROOT}/$ENV{EDEP_TARGET}/lib")
  include_directories("$ENV{EDEP_ROOT}/$ENV{EDEP_TARGET}/include/edep-sim")
  message("EDepSim_DIR is " $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET} )
endif()

#Directory where configuration files like PDFs can be stored.
#add_definitions(-DINSTALL_CONFIG_DIR="${CMAKE_INSTALL_PREFIX}/conf/")

#Write a setup script to set LD_LIBRARY_PATH based on the location of this script.
file(WRITE "setup.sh" "export LD_LIBRARY_PATH=\"`pwd`/`dirname \$BASH_SOURCE`/../lib:\${LD_LIBRARY_PATH}\n\"")
file(APPEND "setup.sh" "export EDEPNEUTRONS_CONF_PATH=\"`pwd`/`dirname \$BASH_SOURCE`/../conf/\"")
#TODO: Handle absolute paths to bash script correctly by not adding pwd.
#file(WRITE  "setup.sh" "BIN_DIR=`dirname \$BASH_SOURCE`\n")
#file(APPEND "setup.sh" "if [[\$BIN_DIR != /*]]\n")
#file(APPEND "setup.sh" "then BIN_DIR=`pwd`/\$BIN_DIR\n")
#file(APPEND "setup.sh" "fi\n")
#file(APPEND "setup.sh" "export LD_LIBRARY_PATH=`pwd`/`dirname $BASH_SOURCE`/../lib\n")
install(FILES "setup.sh" DESTINATION bin)

#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS Core Geom Tree RIO TreePlayer EG)
include(${ROOT_USE_FILE})

#YAML parser library for configuration
find_package( yaml-cpp REQUIRED )

#Make sure that all project files know about all other project files by their 
#path from this directory.
include_directories( "${PROJECT_SOURCE_DIR}" )

#Set up components that used to live in util
add_subdirectory(Base)
add_subdirectory(ROOT)
add_subdirectory(IO)

add_subdirectory(persistency)
add_subdirectory(alg)
add_subdirectory(reco)
add_subdirectory(ana)
add_subdirectory(app)
add_subdirectory(grid)
add_subdirectory(conf)

#Make the results of this build into a package.  Designed to be distributed as a .tar.gz
#Learned to do this from http://agateau.com/2009/cmake-and-make-dist/
set( CPACK_PACKAGE_VERSION_MAJOR "2" )
set( CPACK_PACKAGE_VERSION_MINOR "0" )
set( CPACK_PACKAGE_VERSION_PATCH "0" )
set( CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
set( CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}")
set( CPACK_GENERATOR "TGZ" )

include( CPack )
#Whatever magic this custom target does allows my target pnfs to depend on CPack's files
add_custom_target( myPack COMMAND ${CMAKE_MAKE_PROGRAM} package )

#The dCache file system at Fermilab (/pnfs) does not allow files to be overwritten, but tar files that go to the grid need to be 
#put there.  So, explicitly delete the old tar gile it it exists before installing the new one.  
#
#make pnfs
#
set( FILE_TO_DELETE "${PNFS_DIR}/${CMAKE_PROJECT_NAME}.tar.gz" )
add_custom_target( pnfs COMMAND rm -f ${FILE_TO_DELETE} 
                        COMMAND cp ${CMAKE_PROJECT_NAME}.tar.gz ${FILE_TO_DELETE} 
                        COMMENT "Explicitly removing ${FILE_TO_DELETE} for /pnfs compatibility..." )
add_dependencies( pnfs myPack )
