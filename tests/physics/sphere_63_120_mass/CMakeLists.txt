
add_test(NAME sphere_mass COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_mass.ffea)
set_tests_properties(sphere_mass PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME sphere_mass_check COMMAND ${PYTHON_EXECUTABLE} averages.py)




