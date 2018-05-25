include_directories(${PROJECT_SOURCE_DIR})

#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core)
include(${ROOT_USE_FILE})

add_subdirectory(Style)

#Set up a library to be imported by other projects
file( GLOB Util_ROOT_Base_src *.cxx )
file( GLOB Util_ROOT_Base_headers *.h )

add_library( Util_ROOT_Base ${Util_ROOT_Base_src} )
install( TARGETS Util_ROOT_Base DESTINATION lib EXPORT Util_ROOT_Base_exp )
install( EXPORT Util_ROOT_Base_exp DESTINATION lib/Util_ROOT_Base )
install( FILES ${Util_ROOT_Base_headers} DESTINATION include/ROOT/Base )
target_link_libraries( Util_ROOT_Base ${ROOT_LIBRARIES} Util_Base )
#export( EXPORT Util_ROOT_Base_exp )
