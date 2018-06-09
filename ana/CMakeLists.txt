#Add libraries of plugins
add_library( ana SHARED Analyzer.cpp FSNeutrons.cpp NeutronCand.cpp BirksValidation.cpp NeutronTOF.cpp CandRecoStats.cpp CandTOF.cpp )
target_link_libraries( ana persistency ${EDepSimIO} ${ROOT_LIBRARIES} Util_ROOT_Base Util_Base Truth)
install( TARGETS ana DESTINATION lib )
install( FILES Analyzer.h FSNeutrons.h NeutronCand.h BirksValidation.h NeutronTOF.h CandRecoStats.h CandTOF.h DESTINATION include/ana )
