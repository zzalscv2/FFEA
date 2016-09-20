
add_subdirectory("make_cuboid_mesh")

install(FILES convert_tet_to_net.py cull_small_interior_elements.py
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_volume_tools")

