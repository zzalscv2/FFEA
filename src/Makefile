dir=../bin/
test=../test/
opt = -I../include/ -DFFEA_PARALLEL_WITHIN_BLOB -finline-functions -Wall -Wstrict-aliasing=2 -D_FILE_OFFSET_BITS=64
##-Winline 

GNU_C = g++ 
GNU_CFLAGS = -fopenmp $(opt)
##-march=opteron
INTEL_C = icpc
INTEL_CFLAGS = -openmp -openmp-lib=compat -wd981 -static-intel -static -i-static -limf $(opt)

###############################################
## Set this to either $(GNU_C) or $(INTEL_C) ##
CC = $(GNU_C)
CFLAGS = $(GNU_CFLAGS)
###############################################

OBJECTS = FFEA_return_codes.o BEM_Poisson_Boltzmann.o mesh_node.o \
	SimulationParams.o tetra_element_linear.o \
	SparseSubstitutionSolver.o Face.o Blob.o \
	BiCGSTAB_solver.o CG_solver.o ConjugateGradientSolver.o \
        GaussianQuadrature_1d.o GaussianQuadrature_tri.o \
	MassLumpedSolver.o MassMatrixQuadratic.o MatrixFixedSparsityPattern.o \
	LJ_matrix.o mat_vec_fns.o MersenneTwister.o \
	NearestNeighbourLinkedListCube.o NoMassCGSolver.o \
	PoissonMatrixQuadratic.o SparseMatrixFixedPattern.o \
	SparseMatrixTypes.o SparseMatrixUnknownPattern.o SparsityPattern.o \
	Spring.o World.o VdW_solver.o
SOURCE = FFEA_return_codes.cpp BEM_Poisson_Boltzmann.cpp mesh_node.cpp \
	SimulationParams.cpp tetra_element_linear.cpp \
	SparseSubstitutionSolver.cpp Face.cpp Blob.cpp \
	BiCGSTAB_solver.cpp CG_solver.cpp ConjugateGradientSolver.cpp \
        GaussianQuadrature_1d.cpp GaussianQuadrature_tri.cpp \
	MassLumpedSolver.cpp MassMatrixQuadratic.cpp MatrixFixedSparsityPattern.cpp \
	LJ_matrix.cpp mat_vec_fns.cpp MersenneTwister.cpp \
	NearestNeighbourLinkedListCube.cpp NoMassCGSolver.cpp \
	PoissonMatrixQuadratic.cpp SparseMatrixFixedPattern.cpp \
	SparseMatrixTypes.cpp SparseMatrixUnknownPattern.cpp SparsityPattern.cpp \
	Spring.cpp World.cpp VdW_solver.cpp

all: ffea.o $(OBJECTS) 
	$(CC) $(CFLAGS) -O3 -lm ffea.o $(OBJECTS) -o ffea

debug: ffea.cpp
	$(CC) $(CFLAGS) -ggdb -O0 -lm ffea.cpp -o $(dir)ffea-debug

test: ffea.cpp
	$(CC) $(CFLAGS) -O3 -lm ffea.cpp -o $(test)ffea-test

%.o: %.cpp
	${CC} ${CFLAGS} -c $< -o $@

clean:
	rm -f *.o

install:
	mv ffea $(dir)
