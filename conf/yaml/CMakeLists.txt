file(GLOB defaultConfigs *.yaml)

install(FILES ${defaultConfigs} DESTINATION conf/yaml)
