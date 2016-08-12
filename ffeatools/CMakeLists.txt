
# 1 - Configure and install the ffeatools script:
configure_file("ffeatools.in"
               "${PROJECT_BINARY_DIR}/ffeatools/ffeatools" @ONLY)
install(PROGRAMS "${PROJECT_BINARY_DIR}/ffeatools/ffeatools" 
        DESTINATION  "${CMAKE_INSTALL_PREFIX}/bin")

# 2 - Install python files here:
install(FILES "ffeatools.py"
        DESTINATION ${PYTHONSTUFF})

# 3 - and go through the subdirectories:
add_subdirectory("FFEA_analysis")

add_subdirectory("FFEA_initialise")

# 4 - install the PyMOL plugin:
file (MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/ffeatools/ffea")
set(PYMOL_FILES FFEA_analysis/pymol_plugin/Blob.py
                FFEA_analysis/pymol_plugin/__init__.py
		FFEA_analysis/pymol_plugin/system_names_pokemon.txt
		FFEA_analysis/pymol_plugin/system_names_dbzcharacters.txt
		FFEA_analysis/pymol_plugin/system_names_greekletters.txt
                FFEA_analysis/FFEA_python_modules/FFEA_trajectory.py 
                FFEA_analysis/FFEA_python_modules/FFEA_measurement.py 
                FFEA_analysis/FFEA_python_modules/FFEA_frame.py
                FFEA_analysis/FFEA_python_modules/FFEA_universe.py
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
                FFEA_initialise/FFEA_python_modules/FFEA_pdb.py 
                FFEA_initialise/FFEA_python_modules/FFEA_kinetic_map.py
                FFEA_initialise/FFEA_python_modules/FFEA_kinetic_rates.py
                FFEA_initialise/FFEA_python_modules/FFEA_kinetic_states.py)
file ( COPY ${PYMOL_FILES} DESTINATION "${PROJECT_BINARY_DIR}/ffeatools/ffea" )
set ( FFEA_TAR_BALL "${PROJECT_BINARY_DIR}/ffeatools/ffea.tar.gz")
set ( COMPRESSING TRUE )

# compress the folder.
if (CMAKE_VERSION VERSION_LESS "3.2")
  if (NOT UNIX)
     set ( COMPRESSING FALSE)
     install(DIRECTORY ${PROJECT_BINARY_DIR}/pymol_plugin/ DESTINATION "${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea")
     message(STATUS "PyMOL plugin won't be compressed. The sources will be found in ${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea. You will need to zip or tar.gz the folder yourself before installing it into PyMOL")
  endif(NOT UNIX)
endif (CMAKE_VERSION VERSION_LESS "3.2")

if ( COMPRESSING )
  message(STATUS "PyMOL plugin will be compressed")
  add_custom_target(create_tar ALL COMMAND
                  ${CMAKE_COMMAND} -E tar "cfz" ${FFEA_TAR_BALL} "${PROJECT_BINARY_DIR}/ffeatools/ffea")
  install (FILES ${FFEA_TAR_BALL} DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ffea/plugins/pymol")
endif ( COMPRESSING )

