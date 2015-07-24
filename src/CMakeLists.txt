include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(ffea ffea.cpp
<<<<<<< HEAD
    FFEA_return_codes.cpp BEM_Poisson_Boltzmann.cpp mesh_node.cpp 
=======
    FFEA_return_codes.cpp FFEA_input_reader.cpp 
    BEM_Poisson_Boltzmann.cpp mesh_node.cpp 
>>>>>>> ae21ff1c06f373017017b6c762ab58777e9d28b5
    SimulationParams.cpp tetra_element_linear.cpp 
    SparseSubstitutionSolver.cpp Face.cpp Blob.cpp
    BiCGSTAB_solver.cpp CG_solver.cpp ConjugateGradientSolver.cpp 
    GaussianQuadrature_1d.cpp GaussianQuadrature_tri.cpp 
    MassLumpedSolver.cpp MassMatrixQuadratic.cpp MatrixFixedSparsityPattern.cpp 
<<<<<<< HEAD
    LJ_matrix.cpp mat_vec_fns.cpp MersenneTwister.cpp 
    NearestNeighbourLinkedListCube.cpp NoMassCGSolver.cpp 
    PoissonMatrixQuadratic.cpp SparseMatrixFixedPattern.cpp 
    SparseMatrixTypes.cpp SparseMatrixUnknownPattern.cpp SparsityPattern.cpp 
    Spring.cpp World.cpp VdW_solver.cpp)
=======
    LJ_matrix.cpp BindingSite.cpp mat_vec_fns.cpp MersenneTwister.cpp 
    NearestNeighbourLinkedListCube.cpp NoMassCGSolver.cpp 
    PoissonMatrixQuadratic.cpp SparseMatrixFixedPattern.cpp 
    SparseMatrixTypes.cpp SparseMatrixUnknownPattern.cpp SparsityPattern.cpp 
    Spring.cpp World.cpp VdW_solver.cpp KineticState.cpp KineticBindingSite.cpp)
>>>>>>> ae21ff1c06f373017017b6c762ab58777e9d28b5

target_link_libraries(ffea ${LIBS_TO_LINK})

install (TARGETS ffea RUNTIME DESTINATION bin)

