
add_subdirectory("make_cuboid_mesh")

install(FILES volume_translation.py convert_tet_to_net.py
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Volume_tools")

