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

#ifndef PRECOMP_H_INCLUDED
#define PRECOMP_H_INCLUDED

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include "FFEA_return_codes.h"
#include "dimensions.h"
#include "FFEA_user_info.h"


#include "mat_vec_types.h"
// WARNING: Blob.h will be included after defining PreComp_params! 

using namespace std;
/**
 * @detail 
 * vector<string> types: types of beads present. \n
 * string folder: folder containing the tables. It can be either absolute or relative.\n
 * int inputData: 1 means read .force and .pot files,
 *                 while 2 means read .pot and calculate the forces \n
 */
struct PreComp_params {
  vector<string> types; ///< types of beads present 
  string folder; ///< folder containing the tables. It can be either absolute or relative to the folder containing the ffea input file .
  int inputData; ///< 1 means read .force and .pot files, while 2 means read .pot and calculate the forces
  scalar dist_to_m;
  scalar E_to_J;
};


#include "Blob.h"

class PreComp_solver{
public:
  PreComp_solver();
  ~PreComp_solver();
  int init(PreComp_params *pc_params, SimulationParams *params, Blob **blob_array);
  int solve();
  void reset_fieldenergy(); 
  scalar get_U(scalar x, int typei, int typej);
  scalar get_F(scalar x, int typei, int typej);
  scalar get_field_energy(int index0, int index1);

private: 
  /** msgc and msg are helpful while developing */
  int msgc;
  int msg(string whatever); 
  int msg(int whatever); 
  
  int read_tabulated_values(PreComp_params &pc_params, string kind, scalar *Z, scalar scale_Z);

  int calc_force_from_pot();
  
  scalar finterpolate(scalar *Z, scalar x, int typei, int typej);

  int compute_bead_positions();

  /** delta x in tabulated potentials and forces"  */
  scalar Dx; 
  /** x_range */ 
  scalar x_range[2];
  /** squared x_range */ 
  scalar x_range2[2];
  /** number of pre-computed values per table */ 
  int n_values; 
  /** total number of type interactions */
  int nint; 
  /** number of types of "beads */
  int ntypes; 
  /** pointer to array containing all the values for all the pair potentials. */
  scalar *U;
  /** pointer to array containing all the values for all the pair forces. */
  scalar *F;
  /** interacting elements */
  typedef tetra_element_linear* TELPtr;
  TELPtr *b_elems;
  /** bead types */
  int *b_types; 
  /** number of beads */ 
  int n_beads; 
  /** relative position of the beads to the element they belong, xyzxyzxyz... */
  scalar *b_rel_pos;
  /** absolute position of the beads */
  scalar *b_pos; 
  /** bool "matrix" (array) storing for every pair if it is active or not */
  bool *isPairActive;

  /** Variables to store the energy field data between each pair of blobs */
  scalar **fieldenergy;
  int num_blobs;

};

#endif
