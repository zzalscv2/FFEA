
add_subdirectory("Euler_characteristic_analysis")
add_subdirectory("FFEA_pyPca")
add_subdirectory("FFEA_traj_tools")
add_subdirectory("FFEA_python_modules")

install(DIRECTORY Diffusion_analysis FFEA_gromacs FFEA_meas_tools
                  FFEA_viewer Timestep_analysis
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")

install (FILES FFEA_build_analysis_environment.py FFEA_plot_all_measurements.py 
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")


