
# 1 - Configure and install the FFEA_tools script:
configure_file("FFEA_tools.in"
               "${PROJECT_BINARY_DIR}/FFEA_tools/FFEA_tools" @ONLY)
install(PROGRAMS "${PROJECT_BINARY_DIR}/FFEA_tools/FFEA_tools" 
        DESTINATION  "${CMAKE_INSTALL_PREFIX}/bin")

# 2 - Install python files here:
install(FILES "FFEA_tools.py"
        DESTINATION ${PYTHONSTUFF})

# 3 - and go through the subdirectories:
add_subdirectory("FFEA_analysis")

add_subdirectory("FFEA_initialise")

# 4 - install the PyMOL plugin:
file (MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/ffea")
set(PYMOL_FILES FFEA_analysis/pymol_plugin/Blob.py
                FFEA_analysis/pymol_plugin/__init__.py
                FFEA_analysis/FFEA_python_modules/FFEA_trajectory.py 
                FFEA_analysis/FFEA_python_modules/FFEA_measurement.py 
                FFEA_analysis/FFEA_python_modules/FFEA_frame.py
                FFEA_initialise/FFEA_python_modules/FFEA_binding_sites.py 
                FFEA_initialise/FFEA_python_modules/FFEA_script.py 
                FFEA_initialise/FFEA_python_modules/FFEA_pin.py 
                FFEA_initialise/FFEA_python_modules/FFEA_springs.py 
                FFEA_initialise/FFEA_python_modules/FFEA_material.py 
                FFEA_initialise/FFEA_python_modules/FFEA_topology.py 
                FFEA_initialise/FFEA_python_modules/FFEA_surface.py
                FFEA_initialise/FFEA_python_modules/FFEA_node.py 
                FFEA_initialise/FFEA_python_modules/FFEA_stokes.py 
                FFEA_initialise/FFEA_python_modules/FFEA_vdw.py 
                FFEA_initialise/FFEA_python_modules/FFEA_pdb.py)
file ( COPY ${PYMOL_FILES} DESTINATION "${PROJECT_BINARY_DIR}/ffea" )



