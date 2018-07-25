#include <iostream>
#include <stdlib.h>
#include <omp.h>

#include "FFEA_return_codes.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
//#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"

// ffea_test allows you to build unit tests in FFEA that can be accessed via a .ffeatest script file
// The reason is that otherwise, you'd have to compile a new binary, and you can't easily separate
// out a small section of FFEA codebase to compile - you basically have to compile the whole binary
// each time.

struct ffea_test{

    static int do_ffea_test(std::string filename);

    static int connection_test();

};
