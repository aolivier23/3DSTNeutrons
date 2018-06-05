#Make sure system libraries can be found
link_directories( /usr/local/lib )

add_library(Factory Factory.cpp)
target_link_libraries(Factory)
install(TARGETS Factory DESTINATION lib)

add_executable(EdepApp EdepApp.cpp)
target_link_libraries(EdepApp persistency reco ana ${ROOT_LIBRARIES} yaml-cpp Util_ROOT_Base Util_IO_File ${EDepSimIO})
install(TARGETS EdepApp DESTINATION bin)
