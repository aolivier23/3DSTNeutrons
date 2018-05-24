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
#Force CMake to link paths "correctly".  See the discussion at https://cmake.org/pipermail/cmake/2011-June/044807.html
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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

#Make sure that all project files know about all other project files by their 
#path from this directory.
include_directories( "${PROJECT_SOURCE_DIR}" )

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
