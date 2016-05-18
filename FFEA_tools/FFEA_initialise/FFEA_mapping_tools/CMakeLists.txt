add_executable(make_structure_map make_structure_map.cpp)

install (TARGETS make_structure_map 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")


install(FILES FFEA_map_PDB_to_FFEA.py FFEA_convert_kinetic_map_to_sparse.py FFEA_make_structures_overlap.py FFEA_generate_kinetic_maps.py Kinetic_FFEA_map_apply_to_FFEA_node.py Kinetic_FFEA_map_apply_to_FFEA_map.py Kinetic_FFEA_map_apply_to_PDB.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")
