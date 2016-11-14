#ifndef LJMATRIXHINCLUDED
#define LJMATRIXHINCLUDED

#include <stdio.h>
#include <cstring> 
#include <string>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "dimensions.h"

#define LJI(A,B)	((A) * num_vdw_face_types + (B))

using namespace std;

class LJ_pair {
public:
    scalar vdw_eps;
    scalar vdw_r_eq;
};

class LJ_matrix {
public:
    LJ_matrix();

    int init(string vdw_params_fname, string vdw_type);

    void get_LJ_params(int type1, int type2, scalar *vdw_eps, scalar *vdw_r_eq);

    int get_num_types();

private:
    int init_lj(string vdw_params_fname);
    int init_steric(); 
    LJ_pair *params;
    int num_vdw_face_types;
};


#endif
