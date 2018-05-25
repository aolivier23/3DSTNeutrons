#Make sure ROOT is included
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOT_SYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core)
include(${ROOT_USE_FILE})

add_library(Util_ROOT_Base_Style SelectStyle.cxx DebugStyle.cxx StandardStyle.cxx)
target_link_libraries(Util_ROOT_Base_Style ${ROOT_LIBRARIES})

install(TARGETS Util_ROOT_Base_Style DESTINATION lib EXPORT Util_ROOT_Base_Style_exp)
install( EXPORT Util_ROOT_Base_Style_exp DESTINATION lib/Util_ROOT_Base_Base_Style )
