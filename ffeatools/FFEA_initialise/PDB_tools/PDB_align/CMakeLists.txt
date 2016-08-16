add_executable(align_point_cloud align_point_cloud.cpp)

install(TARGETS align_point_cloud 
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools/PDB_align")

install(FILES align_point_cloud.py rotate_pdb.py
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/PDB_tools/PDB_align")
