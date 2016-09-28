add_executable(build_euler_characteristic_table build_euler_characteristic_table.cpp)

install(TARGETS build_euler_characteristic_table
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/Euler_characteristic_analysis")
