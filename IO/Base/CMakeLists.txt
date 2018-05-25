#Set path where executables will be made
#set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin) 

include_directories(${PROJECT_SOURCE_DIR})

link_directories( /usr/local/lib )

add_library( Util_IO_Base decode.cxx )
install( TARGETS Util_IO_Base DESTINATION lib EXPORT Util_IO_Base )
target_link_libraries( Util_IO_Base Util_Base )
install( EXPORT Util_IO_Base DESTINATION lib/Util_IO_Base )
install( FILES decode.cxx DESTINATION include/IO/Base )
#export( EXPORT Util_IO_Base ) #TODO: convert this to CMake 2 so that I can use it on the minervagpvms
