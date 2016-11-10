
add_subdirectory("Euler_characteristic_analysis")
add_subdirectory("FFEA_pyPca")
add_subdirectory("FFEA_traj_tools")

install (FILES FFEA_build_analysis_environment.py __init__.py
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")


