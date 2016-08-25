
# 1 - Configure and install the FFEA_tools script:
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

add_subdirectory("modules")

# 4 - install the PyMOL plugin:
set ( FFEA_TAR_DIR "${PROJECT_BINARY_DIR}/ffeatools")
file (MAKE_DIRECTORY "${FFEA_TAR_DIR}/FFEAplugin")
set(PYMOL_FILES FFEA_analysis/pymol_plugin/Blob.py
                FFEA_analysis/pymol_plugin/__init__.py
		FFEA_analysis/pymol_plugin/system_names_pokemon.txt
		FFEA_analysis/pymol_plugin/system_names_dbzcharacters.txt
		FFEA_analysis/pymol_plugin/system_names_greekletters.txt
                modules/FFEA_trajectory.py 
                modules/FFEA_measurement.py 
                modules/FFEA_frame.py
                modules/FFEA_universe.py
                modules/FFEA_binding_sites.py 
                modules/FFEA_script.py 
                modules/FFEA_pin.py 
                modules/FFEA_springs.py 
                modules/FFEA_material.py 
                modules/FFEA_topology.py 
                modules/FFEA_surface.py
                modules/FFEA_node.py 
                modules/FFEA_stokes.py 
                modules/FFEA_vdw.py 
                modules/FFEA_lj.py 
                modules/FFEA_pdb.py 
                modules/FFEA_kinetic_map.py
                modules/FFEA_kinetic_rates.py
                modules/FFEA_kinetic_states.py)
file ( COPY ${PYMOL_FILES} DESTINATION "${FFEA_TAR_DIR}/FFEAplugin" )
set ( FFEA_TAR_BALL "${FFEA_TAR_DIR}/FFEAplugin.tar.gz")
set ( COMPRESSING TRUE )

# compress the folder.
if (CMAKE_VERSION VERSION_LESS "3.2")
  if (NOT UNIX)
     set ( COMPRESSING FALSE)
     install(DIRECTORY ${FFEA_TAR_DIR}/FFEAplugin/ DESTINATION "${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea")
     message(STATUS "PyMOL plugin won't be compressed. The sources will be found in ${CMAKE_INSTALL_PREFIX}/share/plugins/pymol/ffea. You will need to zip or tar.gz the folder yourself before installing it into PyMOL")
  endif(NOT UNIX)
endif (CMAKE_VERSION VERSION_LESS "3.2")

if ( COMPRESSING )
  message(STATUS "PyMOL plugin will be compressed")
  add_custom_target(create_tar ALL COMMAND
                  ${CMAKE_COMMAND} -E tar "cfz" ${FFEA_TAR_BALL} "${FFEA_TAR_DIR}/FFEAplugin/")
  install (FILES ${FFEA_TAR_BALL} DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ffea/plugins/pymol")
endif ( COMPRESSING )

