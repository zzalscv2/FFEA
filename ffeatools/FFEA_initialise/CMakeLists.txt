
add_subdirectory("EM_density_map_tools")
add_subdirectory("FFEA_mapping_tools")
add_subdirectory("Surface_tools")
add_subdirectory("FFEA_volume_tools")

install(DIRECTORY FFEA_convert_from_volume
        FFEA_material_tools FFEA_pdb_tools FFEA_pin_tools
        Kinetic_FFEA_convert_from_FFEA FFEA_vdw_tools
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise")

install (FILES FFEA_build_simulation_environment.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise")

