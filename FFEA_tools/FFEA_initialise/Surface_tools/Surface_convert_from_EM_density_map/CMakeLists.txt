add_executable(emdb_map_to_ffea emdb_map_to_ffea.cpp)

install(TARGETS emdb_map_to_ffea
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map")
