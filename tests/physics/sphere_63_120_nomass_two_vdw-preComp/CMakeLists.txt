
add_test(NAME precomp COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_two.ffea)
set_tests_properties(precomp PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME precomp_check COMMAND ${PYTHON_EXECUTABLE} check.py)
set_tests_properties(precomp_check PROPERTIES DEPENDS precomp)




