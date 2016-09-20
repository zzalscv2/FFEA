
add_subdirectory("Euler_characteristic_analysis")
add_subdirectory("FFEA_pyPca")
add_subdirectory("FFEA_traj_tools")

install(DIRECTORY FFEA_meas_tools FFEA_viewer
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")

install (FILES FFEA_build_analysis_environment.py 
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")


