dir=../bin/
test=../test/
opt = -finline-functions -Wall -Wstrict-aliasing=2 -D_FILE_OFFSET_BITS=64
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
	SparseSubstitutionSolver.o Face.o Blob.o
SOURCE = FFEA_return_codes.cpp BEM_Poisson_Boltzmann.cpp mesh_node.cpp \
	SimulationParams.cpp tetra_element_linear.cpp \
	SparseSubstitutionSolver.cpp Face.cpp Blob.cpp

all: ffea.o $(OBJECTS) 
	$(CC) $(CFLAGS) -O3 -lm ffea.o $(OBJECTS) -o $(dir)ffea

debug: ffea.cpp
	$(CC) $(CFLAGS) -ggdb -O0 -lm ffea.cpp -o $(dir)ffea-debug

test: ffea.cpp
	$(CC) $(CFLAGS) -O3 -lm ffea.cpp -o $(test)ffea-test

%.o: %.cpp
	${CC} ${CFLAGS} -c $< -o $@

clean:
	rm -f *.o
