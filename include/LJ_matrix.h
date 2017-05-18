// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

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
    LJ_pair();
    ~LJ_pair();
    scalar vdw_eps;
    scalar vdw_r_eq;
};

class LJ_matrix {
public:
    LJ_matrix();
    ~LJ_matrix(); 

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
