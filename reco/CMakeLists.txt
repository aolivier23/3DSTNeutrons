#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Get util library
message( "Including util library from ${UTIL_DIR}\n" )
include_directories( "${UTIL_DIR}/include" )
link_directories("${UTIL_DIR}/lib" )

add_subdirectory(alg)

add_library( reco SHARED Reconstructor.cpp NeutronHits.cpp NoGridNeutronHits.cpp GridNeutronHits.cpp 
                         AdjacentClusters.cpp MergedClusters.cpp TreeNeutronHits.cpp GridAllHits.cpp 
                         CCQEChargedFSFilter.cpp)
target_link_libraries( reco persistency ${ROOT_LIBRARIES} Util_ROOT_Base Util_Base Truth RecoAlgs Geo)
install( TARGETS reco DESTINATION lib )
install( FILES Reconstructor.h NeutronHits.h NoGridNeutronHits.h GridNeutronHits.h AdjacentClusters.h 
               MergedClusters.h TreeNeutronHits.h GridAllHits.h CCQEChargedFSFilter.h 
        DESTINATION include/reco )
