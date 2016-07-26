
add_test(NAME sphere_nomass COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass.ffea)
set_tests_properties(sphere_nomass PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME sphere_nomass_check COMMAND ${PYTHON_EXECUTABLE} averages.py)
set_tests_properties(sphere_nomass_check PROPERTIES DEPENDS sphere_nomass)




