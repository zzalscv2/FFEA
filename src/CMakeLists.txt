# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(ffea ffea.cpp
    FFEA_return_codes.cpp FFEA_input_reader.cpp 
    BEM_Poisson_Boltzmann.cpp mesh_node.cpp 
    SimulationParams.cpp tetra_element_linear.cpp 
    SparseSubstitutionSolver.cpp Face.cpp Blob.cpp
    BiCGSTAB_solver.cpp CG_solver.cpp ConjugateGradientSolver.cpp 
    GaussianQuadrature_1d.cpp GaussianQuadrature_tri.cpp 
    MassLumpedSolver.cpp MassMatrixLinear.cpp MassMatrixQuadratic.cpp MatrixFixedSparsityPattern.cpp 
    LJ_matrix.cpp BindingSite.cpp mat_vec_fns.cpp mat_vec_fns_II.cpp
    NearestNeighbourLinkedListCube.cpp NoMassCGSolver.cpp 
    PoissonMatrixQuadratic.cpp SparseMatrixFixedPattern.cpp 
    SparseMatrixTypes.cpp SparseMatrixUnknownPattern.cpp SparsityPattern.cpp 
    Spring.cpp World.cpp VdW_solver.cpp 
    Steric_solver.cpp LJSteric_solver.cpp KineticState.cpp PreComp_solver.cpp 
    CheckTetrahedraOverlap_II.cpp VolumeIntersection.cpp
    FFEA_user_info.cpp RngStream.cpp )
    
if (USE_BOOST_INTERNAL)
   add_dependencies(ffea ourboost)
endif(USE_BOOST_INTERNAL)
if (USE_EIGEN3_INTERNAL)
   add_dependencies(ffea getEigen)
endif(USE_EIGEN3_INTERNAL)


target_link_libraries(ffea ${LIBS_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})

install (TARGETS ffea RUNTIME DESTINATION bin)

