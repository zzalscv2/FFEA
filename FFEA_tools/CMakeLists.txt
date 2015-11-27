
# 1 - Configure and install the FFEA_tools script:
configure_file("FFEA_tools.in"
               "${PROJECT_BINARY_DIR}/FFEA_tools/FFEA_tools" @ONLY)
install(PROGRAMS "${PROJECT_BINARY_DIR}/FFEA_tools/FFEA_tools" 
        DESTINATION  "${CMAKE_INSTALL_PREFIX}/bin")

# 2 - Install python files here:
install(FILES "FFEA_tools.py"
        DESTINATION ${PYTHONSTUFF})

# 3 - and go through the subdirectories:
add_subdirectory("FFEA_initialise")

add_subdirectory("FFEA_analysis")
