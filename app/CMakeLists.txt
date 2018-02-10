#Make sure system libraries can be found
link_directories( /usr/local/lib )

#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Get util library
#TODO: Make util a real package so I don't need to do this
include_directories( "${UTIL_DIR}/src" )
link_directories("${UTIL_DIR}/lib" )

add_library(Factory Factory.cpp)
target_link_libraries(Factory)
install(TARGETS Factory DESTINATION lib)

add_executable(EdepApp EdepApp.cpp)
target_link_libraries(EdepApp persistency reco ana ${ROOT_LIBRARIES} Util_IO_Option_runtime_CmdLine Util_IO_Option_runtime_Policy 
                      Util_IO_Option_runtime_Options Util_Base Util_IO_Option_runtime_PolicyBase 
                      Util_IO_Option_PrintCmdLine Util_ROOT_Base Util_IO_File edepsim)
install(TARGETS EdepApp DESTINATION bin)
