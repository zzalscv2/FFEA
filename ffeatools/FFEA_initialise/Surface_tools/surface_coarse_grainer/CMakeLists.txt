add_executable(surface_coarse_grainer_final surface_coarse_grainer_final.cpp)

install(TARGETS surface_coarse_grainer_final
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Surface_tools/surface_coarse_grainer")

install(FILES __init__.py surface_coarse_grainer_final.py 
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Surface_tools/surface_coarse_grainer")
