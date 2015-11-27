add_executable(make_cuboid_mesh make_cuboid_mesh.c)

install(TARGETS make_cuboid_mesh
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Volume_tools/make_cuboid_mesh")
