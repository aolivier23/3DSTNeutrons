link_directories( /usr/local/lib )

#Set up a library to be imported by other projects
#file( GLOB Util_IO_Option_runtime_src *.cxx )

file( GLOB Util_IO_Option_runtime_headers *.h )

include_directories( "${PROJECT_SOURCE_DIR}" )
#add_library( Util_IO_Option_runtime ${Util_IO_Option_runtime_src} )
add_library( Util_IO_Option_runtime_PolicyBase Policy.cxx )
add_library( Util_IO_Option_runtime_Policy ExactlyOnce.cxx Accumulate.cxx Counter.cxx Exists.cxx Help.cxx )
add_library( Util_IO_Option_runtime_Options Options.cxx )
add_library( Util_IO_Option_runtime_CmdLine CmdLine.cxx )

install( TARGETS Util_IO_Option_runtime_PolicyBase DESTINATION lib EXPORT Util_IO_Option_runtime_PolicyBase_exp )
install( TARGETS Util_IO_Option_runtime_Policy DESTINATION lib EXPORT Util_IO_Option_runtime_Policy_exp )
install( TARGETS Util_IO_Option_runtime_Options DESTINATION lib EXPORT Util_IO_Option_runtime_Options_exp )
install( TARGETS Util_IO_Option_runtime_CmdLine DESTINATION lib EXPORT Util_IO_Option_runtime_CmdLine_exp )

target_link_libraries(Util_IO_Option_runtime_PolicyBase Util_Base)
target_link_libraries(Util_IO_Option_runtime_Policy Util_IO_Option_runtime_PolicyBase)
target_link_libraries(Util_IO_Option_runtime_Options Util_IO_Base)
target_link_libraries(Util_IO_Option_runtime_CmdLine Util_IO_Option_runtime_Options Util_IO_Option_runtime_Policy Util_Base )

install( EXPORT Util_IO_Option_runtime_PolicyBase_exp DESTINATION lib/Util_IO_Option_runtime_PolicyBase )
install( EXPORT Util_IO_Option_runtime_CmdLine_exp DESTINATION lib/Util_IO_Option_runtime_CmdLine )
install( EXPORT Util_IO_Option_runtime_Policy_exp DESTINATION lib/Util_IO_Option_runtime_Policy )
install( EXPORT Util_IO_Option_runtime_Options_exp DESTINATION lib/Util_IO_Option_runtime_Options )
install( FILES ${Util_IO_Option_runtime_headers} DESTINATION include/IO/Option/runtime )
#export( EXPORT Util_IO_Option_runtime_exp )

add_executable( RuntimeOptExample RuntimeOptExample.cpp )
target_link_libraries( RuntimeOptExample Util_IO_Option_runtime_CmdLine Util_IO_Option_runtime_Policy Util_IO_Option_runtime_Options )
install( TARGETS RuntimeOptExample DESTINATION bin ) #add stringTest to work with compile-time strings
