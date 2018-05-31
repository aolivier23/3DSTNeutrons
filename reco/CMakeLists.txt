add_subdirectory(alg)

add_library( reco SHARED Reconstructor.cpp NeutronHits.cpp NoGridNeutronHits.cpp GridNeutronHits.cpp 
                         AdjacentClusters.cpp MergedClusters.cpp TreeNeutronHits.cpp GridAllHits.cpp 
                         CCQEChargedFSFilter.cpp CandFromTOF.cpp CandFromPDF.cpp CandFromCluster.cpp)
target_link_libraries( reco persistency ${ROOT_LIBRARIES} Util_ROOT_Base Util_Base Truth RecoAlgs Geo)
install( TARGETS reco DESTINATION lib )
install( FILES Reconstructor.h NeutronHits.h NoGridNeutronHits.h GridNeutronHits.h AdjacentClusters.h 
               MergedClusters.h TreeNeutronHits.h GridAllHits.h CCQEChargedFSFilter.h CandFromTOF.h
               CandFromPDF.h CandFromCluster.h
         DESTINATION include/reco )
