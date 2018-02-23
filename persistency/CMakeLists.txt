# The following CMake code taken directly from Clark McGrew's edepsim since 
# he seems to have done an awesome job.  
# Make sure the current directories are available for the root
# dictionary generation.
include_directories(${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Build the dictionary for the i/o classes.
ROOT_GENERATE_DICTIONARY(G__persistency
                         MCHit.h MCCluster.h
                         OPTIONS -inlineInputHeader
                         LINKDEF LinkDef.h)

add_library(persistency SHARED MCHit.cpp MCCluster.cpp G__persistency.cxx)
target_link_libraries(persistency ${ROOT_LIBRARIES})

install(TARGETS persistency LIBRARY DESTINATION lib)
install(FILES MCHit.h MCCluster.h DESTINATION include/persistency) 

# If this is ROOT6 or later, then install the rootmap and pcm files.
if(${ROOT_VERSION} VERSION_GREATER 6)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libpersistency.rootmap
    DESTINATION lib)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libpersistency_rdict.pcm
    DESTINATION lib)
endif(${ROOT_VERSION} VERSION_GREATER 6)
