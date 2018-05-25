#Set path where executables will be made
#set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin) 

#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR})

link_directories( /usr/local/lib )
link_directories( ${CMAKE_SOURCE_DIR}/Util/IO )

#Set up a library to be imported by other projects
file( GLOB Util_IO_File_src *.cxx )
file( GLOB Util_IO_File_headers *.h )

add_library( Util_IO_File ${Util_IO_File_src} )
install( TARGETS Util_IO_File DESTINATION lib EXPORT Util_IO_File_exp )
target_link_libraries( Util_IO_File Util_IO_Base )
install( EXPORT Util_IO_File_exp DESTINATION lib/Util_IO_File_exp )
install( FILES ${Util_IO_File_headers} DESTINATION include/IO/File )
#export( EXPORT Util_IO_File_exp )
