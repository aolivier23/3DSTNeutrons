#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Get util library
#TODO: Make util a real package so I don't need to do this
set( UTIL_DIR "/dune/app/users/aolivier/nd_neutrons/util")
include_directories( "${UTIL_DIR}/src" )
link_directories( "${UTIL_DIR}/lib" )

#Add libraries of plugins
add_library( ana SHARED Analyzer.cpp FSNeutrons.cpp NeutronCand.cpp )
target_link_libraries( ana persistency edepsim ${ROOT_LIBRARIES} Util_ROOT_Base)
install( TARGETS ana DESTINATION lib )
install( FILES Analyzer.h FSNeutrons.h NeutronCand.h DESTINATION include/ana )
