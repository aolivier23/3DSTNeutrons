link_directories( /usr/local/lib )
include_directories( "${PROJECT_SOURCE_DIR}" )

add_library(Util_IO_Option_PrintCmdLine PrintCmdLine.cxx)
install( TARGETS Util_IO_Option_PrintCmdLine DESTINATION lib )

add_subdirectory( runtime )
