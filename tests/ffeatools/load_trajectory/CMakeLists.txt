set (TESTPYTHONTRAJ "${PROJECT_BINARY_DIR}/tests/ffeatools/load_trajectory")
file (COPY unit_test_traj.ftj DESTINATION ${TESTPYTHONTRAJ})
file (COPY python_load_trajectory.py DESTINATION ${TESTPYTHONTRAJ})
add_test(NAME python_load_trajectory COMMAND ${PYTHON_EXECUTABLE} python_load_trajectory.py)
