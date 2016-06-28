
add_subdirectory("EM_density_map_tools")
add_subdirectory("FFEA_mapping_tools")
add_subdirectory("FFEA_node_tools")
add_subdirectory("PDB_tools")
add_subdirectory("Surface_tools")
add_subdirectory("Volume_tools")

install(DIRECTORY FFEA_convert_from_volume
        FFEA_material_tools FFEA_pdb_tools FFEA_pin_tools
        FFEA_python_modules FFEA_script_tools FFEA_vdw_tools
        Kinetic_FFEA_convert_from_FFEA FFEA_binding_site_tools
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise")

install (FILES FFEA_build_simulation_environment.py FFEA_rerun.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise")

install (FILES FFEA_python_modules/FFEA_binding_sites.py FFEA_python_modules/FFEA_script.py FFEA_python_modules/FFEA_pin.py FFEA_python_modules/FFEA_springs.py FFEA_python_modules/FFEA_material.py FFEA_python_modules/FFEA_topology.py FFEA_python_modules/FFEA_surface.py DESTINATION "${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea")
install (FILES FFEA_python_modules/FFEA_node.py FFEA_python_modules/FFEA_stokes.py FFEA_python_modules/FFEA_vdw.py FFEA_python_modules/FFEA_pdb.py DESTINATION "${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea")
