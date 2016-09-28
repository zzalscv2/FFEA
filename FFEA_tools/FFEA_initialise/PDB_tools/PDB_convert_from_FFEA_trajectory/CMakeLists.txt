add_executable(FFEA_convert_traj_to_pdb FFEA_convert_traj_to_pdb.cpp)

install(TARGETS FFEA_convert_traj_to_pdb
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory")
install(FILES FFEA_convert_traj_to_pdb.py
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory")
