
add_subdirectory("PDB_align")
add_subdirectory("PDB_convert_from_FFEA_trajectory")

install(FILES "PDB_combine.py" "PDB_renumber.py" "PDB_rotate.py"
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools")

install(DIRECTORY "PDB_slice" "PDB_convert_to_FFEA_trajectory"
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools")
