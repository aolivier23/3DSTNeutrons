#Set up a library to be imported by other projects
add_library( Util_Base exception.cxx )
target_link_libraries(Util_Base)

install( TARGETS Util_Base DESTINATION lib)
install( FILES exception.h DESTINATION include/Base )
