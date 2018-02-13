#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

add_library(Truth TruthFunc.cpp)
target_link_libraries(Truth edepsim)
install(TARGETS Truth DESTINATION lib)
install(FILES TruthFunc.h DESTINATION include)
