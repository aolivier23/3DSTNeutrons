add_library(Truth TruthFunc.cpp)
target_link_libraries(Truth ${EDepSimIO})
install(TARGETS Truth DESTINATION lib)
install(FILES TruthFunc.h DESTINATION include)
