add_executable(make_structure_map2 make_structure_map.cpp)

set_target_properties(make_structure_map2 
         PROPERTIES OUTPUT_NAME make_structure_map)

install (TARGETS make_structure_map2 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_node_tools")

install (FILES calc_centroid.py FFEA_node_shear_object.py
	 Generate_rotation_matrix.py node_translation.py transform_node.py
         translate_node_to_node_centroid.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_node_tools")
