add_library(Geo SHARED GeoFunc.cpp)
target_link_libraries(Geo ${ROOT_LIBRARIES})
install(TARGETS Geo DESTINATION lib)

add_library(RecoAlgs SHARED GridHits.cpp Octree.cpp)
target_link_libraries(RecoAlgs Geo ${ROOT_LIBRARIES} ${EDepSimIO})
install(TARGETS RecoAlgs DESTINATION lib)

install(FILES GeoFunc.h GridHits.h DESTINATION include)
