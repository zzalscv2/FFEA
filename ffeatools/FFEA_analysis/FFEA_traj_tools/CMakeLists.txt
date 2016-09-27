add_executable(FFEA_extract_positions FFEA_extract_positions.c)
add_executable(FFEA_extract_specific_nodes FFEA_extract_specific_nodes.c)

install(TARGETS FFEA_extract_positions FFEA_extract_specific_nodes
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

# install(TARGETS FFEA_extract_positions FFEA_extract_specific_nodes
#         DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

install(FILES 
        convert_old_traj_to_new.py Dx_convert_from_xyz_position_trajectory.py
        Dx_smooth_convert_from_xyz_position_trajectory.py FFEA_convert_old_traj_to_new.py
        FFEA_map_trajectory_to_PDB.py FFEA_thin_trajectory.py FFEA_traj_to_nodes.py
        FFEA_traj_to_PDB_traj.py FFEA_convert_traj_to_pdb.py FFEA_trim_trajectory.py FFEA_split_trajectory.py
        FFEA_get_snapshots_in_nodes.py FFEA_strip_equilibration.py FFEA_get_num_frames.py
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

