#Set up to build against ROOT
#I learned how to do this here: https://root.cern.ch/how/integrate-root-my-project-cmake
#list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
#find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core)
#include(${ROOT_USE_FILE})

add_subdirectory(Style)

add_library( Util_ROOT_Base TFileSentry.cxx TCollectionSTLIter.cxx )
target_link_libraries( Util_ROOT_Base ${ROOT_LIBRARIES} Util_Base )
install( TARGETS Util_ROOT_Base DESTINATION lib )
install( FILES TFileSentry.h TCollectionSTLIter.h DESTINATION include/ROOT/Base )
