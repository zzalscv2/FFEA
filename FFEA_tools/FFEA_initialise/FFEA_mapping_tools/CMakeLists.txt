add_executable(make_structure_map make_structure_map.cpp)

install (TARGETS make_structure_map 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")


install(FILES Kinetic_FFEA_map_convert_to_sparse.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")
