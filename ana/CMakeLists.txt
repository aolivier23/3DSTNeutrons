#Set up edep-sim
set(EDEP_DIR $ENV{EDEP_ROOT}/$ENV{EDEP_TARGET})
include_directories( "${EDEP_DIR}/include/edep-sim" )
link_directories( "${EDEP_DIR}/lib" )

#Add libraries of plugins
add_library( ana SHARED Analyzer.cpp FSNeutrons.cpp NeutronCand.cpp BirksValidation.cpp NeutronTOF.cpp CandRecoStats.cpp CandTOF.cpp )
target_link_libraries( ana persistency edepsim ${ROOT_LIBRARIES} Util_ROOT_Base Truth)
install( TARGETS ana DESTINATION lib )
install( FILES Analyzer.h FSNeutrons.h NeutronCand.h BirksValidation.h NeutronTOF.h CandRecoStats.h CandTOF.h DESTINATION include/ana )
