dir=../bin/
test=../test/
opt = -finline-functions -Wall -Wstrict-aliasing=2 -D_FILE_OFFSET_BITS=64
##-Winline 

GNU_C = g++ -fopenmp $(opt)
##-march=opteron
INTEL_C = icpc -openmp -openmp-lib=compat -wd981 -static-intel -static -i-static -limf $(opt)

###############################################
## Set this to either $(GNU_C) or $(INTEL_C) ##
CC = $(GNU_C)
###############################################

all: ffea.cpp
	$(CC) -O3 -lm ffea.cpp -o $(dir)ffea

debug: ffea.cpp
	$(CC) -ggdb -O0 -lm ffea.cpp -o $(dir)ffea

test: ffea.cpp
	$(CC) -O3 -lm ffea.cpp -o $(test)ffea-test

clean:
	rm -f *.o
