add_executable(FFEA_convert_traj_to_mdcrd FFEA_convert_traj_to_mdcrd.cpp)

install(TARGETS FFEA_convert_traj_to_mdcrd
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

install(TARGETS FFEA_convert_traj_to_mdcrd
        DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

install(FILES 
        FFEA_compare_modes.py FFEA_convert_FFEA_to_atomic.py FFEA_get_PCA_animations.py
        FFEA_get_PCA_eigensystem.py FFEA_test_convergence.py plot_eigensystem_comparison.py
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

