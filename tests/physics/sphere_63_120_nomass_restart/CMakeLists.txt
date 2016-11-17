
add_test(NAME restart_run_I COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_10steps.ffea)
set_tests_properties(restart_run_I PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

add_test(NAME restart_run_II COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_1-5steps.ffea)
set_tests_properties(restart_run_II PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

add_test(NAME restart_run_III COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_6-10steps.ffea)
set_tests_properties(restart_run_III PROPERTIES DEPENDS restart_run_II ENVIRONMENT OMP_NUM_THREADS=1)

add_test(NAME restart_check COMMAND ${PYTHON_EXECUTABLE} checkRuns.py)
set_tests_properties(restart_check PROPERTIES DEPENDS "restart_run_I;restart_run_II;restart_run_III")




