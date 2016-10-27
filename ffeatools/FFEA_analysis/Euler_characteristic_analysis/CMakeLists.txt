add_executable(build_euler_characteristic_table build_euler_characteristic_table.cpp)

install(TARGETS build_euler_characteristic_table
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/Euler_characteristic_analysis")

install(FILES __init__.py DESTINATION  "${PYTHONSTUFF}/FFEA_analysis/Euler_characteristic_analysis")

# install(TARGETS build_euler_characteristic_table
#         DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
