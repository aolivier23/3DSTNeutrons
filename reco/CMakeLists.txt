#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Get util library
include_directories( "${UTIL_DIR}/src" )
link_directories("${UTIL_DIR}/lib" )

#Add libraries of plugins
add_library( octree SHARED Octree )
target_link_libraries( octree ${ROOT_LIBRARIES} )
install( FILES Octree.cpp DESTINATION lib )

add_library( reco SHARED Reconstructor.cpp NeutronHits.cpp NoGridNeutronHits.cpp AdjacentClusters.cpp MergedClusters.cpp GridNeutronHits.cpp)
target_link_libraries( reco octree persistency ${ROOT_LIBRARIES} Util_ROOT_Base Util_Base)
install( TARGETS reco DESTINATION lib )
install( FILES Reconstructor.h NeutronHits.h NoGridNeutronHits.h AdjacentClusters.h MergedClusters.h GridNeutronHits.h DESTINATION include/reco )
