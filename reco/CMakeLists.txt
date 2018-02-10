#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Get util library
include_directories( "${UTIL_DIR}/src" )
link_directories("${UTIL_DIR}/lib" )

#Add libraries of plugins
add_library( reco SHARED Reconstructor.cpp NeutronHits.cpp AdjacentClusters.cpp )
target_link_libraries( reco persistency ${ROOT_LIBRARIES} Util_ROOT_Base Util_Base)
install( TARGETS reco DESTINATION lib )
install( FILES Reconstructor.h NeutronHits.h AdjacentClusters.h DESTINATION include/reco )
