# CXXFLAGS="-O0 -ggdb -std=c++11"
CXXFLAGS="-O3 -std=c++11"
CPPFLAGS="-I../../../include"
g++ $CXXFLAGS $CPPFLAGS ../../../src/mat_vec_fns_II.cpp -c -o mvf.o 
g++ $CXXFLAGS $CPPFLAGS ../../../src/BlobLite.cpp -c -o bl.o
g++ $CXXFLAGS $CPPFLAGS ../../../src/CheckTetrahedraOverlap.cpp -c -o cto.o 
g++ $CXXFLAGS $CPPFLAGS ../../../src/VolumeIntersection.cpp -c -o vi.o
g++ $CXXFLAGS $CPPFLAGS ../../../src/FFEA_return_codes.cpp -c -o frc.o
g++ $CXXFLAGS $CPPFLAGS checkIntersections.cpp  -c -o myway.o
g++ $CXXFLAGS frc.o myway.o mvf.o cto.o bl.o vi.o -o a.out
