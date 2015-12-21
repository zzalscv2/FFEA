add_executable(make_structure_map make_structure_map.cpp)

install (TARGETS make_structure_map 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")


install(FILES FFEA_convert_kinetic_map_to_sparse.py FFEA_make_structures_overlap.py FFEA_generate_kinetic_maps.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")
