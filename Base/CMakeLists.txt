#Set path where executables will be made
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin) 

link_directories( /usr/local/lib )

#Set up a library to be imported by other projects
#file( GLOB Util_IO_Option_runtime_src *.cxx )

include_directories( "${PROJECT_SOURCE_DIR}" )
add_library( Util_Base exception.cxx )

install( TARGETS Util_Base DESTINATION lib EXPORT Util_Base_exp )

install( EXPORT Util_Base_exp DESTINATION lib/Util_Base )
install( FILES exception.h DESTINATION include/Base )
