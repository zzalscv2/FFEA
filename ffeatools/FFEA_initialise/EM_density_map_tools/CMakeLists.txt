add_executable(add_maps add_maps.cpp)
add_executable(convert_pdb_to_emdb_map convert_pdb_to_emdb_map.cpp)

install (TARGETS add_maps convert_pdb_to_emdb_map RUNTIME 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/EM_density_map_tools")


