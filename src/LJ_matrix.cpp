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

#include "LJ_matrix.h"
LJ_pair::LJ_pair() {
    Emin = 0;
    Rmin = 0;
}


LJ_pair::~LJ_pair() {
    Emin = 0;
    Rmin = 0;
}

LJ_matrix::LJ_matrix() {
    params = NULL;
    num_ssint_face_types = 0;
}


LJ_matrix::~LJ_matrix() {
    delete[] params; 
    params = NULL;
    num_ssint_face_types = 0;
}


int LJ_matrix::init(string ssint_params_fname, string ssint_type, int calc_ssint) {

    int err;
    // In that case we do not need an input parameter file, 
    //      but just to initialise a number of values
    if (calc_ssint == 0) {
       err = init_steric(); 
    } else {
       err = init_lj(ssint_params_fname); 
    } 
    return err; 
} 

int LJ_matrix::init_steric() {

    num_ssint_face_types = 1;
    // Allocate the memory for that LJ pair
    params = new(std::nothrow) LJ_pair[num_ssint_face_types * num_ssint_face_types];
    if (params == NULL) {
        FFEA_ERROR_MESSG("Unable to allocate memory for LJ matrix.\n")
    }
    params[LJI(0, 0)].Emin = 0; 
    params[LJI(0, 0)].Rmin = 0; 

    return FFEA_OK; 
}

int LJ_matrix::init_lj(string ssint_params_fname) {
    FILE *in = NULL;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(ssint_params_fname.c_str(), "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(ssint_params_fname.c_str())
    }
    printf("\tReading in ssint forcefield parameters file: %s\n", ssint_params_fname.c_str());

    // first line should be the file type "ffea ssint forcefield params file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of ffea ssint forcefield params file\n")
    }
    if (strcmp(line, "ffea vdw forcefield params file\n") != 0 && strcmp(line, "ffea ssint forcefield params file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea ssint forcefield params file' (read '%s') \n", line)
    }

    // read in the number of nodes in the file
    if (fscanf(in, "num_vdw_face_types %d\n", &num_ssint_face_types) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of ssint face types\n")
    }
    printf("\t\tNumber of ssint face types = %d\n", num_ssint_face_types);

    // Allocate the memory for all these LJ pairs
    params = new(std::nothrow) LJ_pair[num_ssint_face_types * num_ssint_face_types];
    if (params == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Unable to allocate memory for LJ matrix.\n")
    }

    // Fill the matrix
    double Emin = 0.0, Rmin = 0.0;
    for (int i = 0; i < num_ssint_face_types; i++) {
        for (int j = 0; j < num_ssint_face_types; j++) {
            if (fscanf(in, " (%le,%le) ", &Emin, &Rmin) != 2) {
                fclose(in);
                FFEA_ERROR_MESSG("Error reading from ffea ssint forcefield params file at (%d,%d)\n", i, j);
            }

            if (Emin <= 0) {
                FFEA_ERROR_MESSG("Required: 'Emin' must be greater than 0 if you wish to use ssint (calc_ssint is 1)\n")
            }
            if (Rmin <= 0) {
                FFEA_ERROR_MESSG("Required: 'Rmin' must be greater than 0 if you wish to use ssint (calc_ssint is 1)\n")
            }

            params[LJI(i, j)].Emin = Emin * mesoDimensions::area * mesoDimensions::area / mesoDimensions::Energy ;
            params[LJI(i, j)].Rmin = Rmin / mesoDimensions::length ;
        }
    }

    fclose(in);

    printf("\t\tRead %d ssint forcefield parameter entries from %s\n", num_ssint_face_types * num_ssint_face_types, ssint_params_fname.c_str());

    return FFEA_OK;
}

void LJ_matrix::get_LJ_params(int type1, int type2, scalar *Emin, scalar *Rmin) {
    if (type1 < 0 || type1 > num_ssint_face_types - 1) {
        printf("Frog1 %d %d\n", type1, num_ssint_face_types - 1);
        return;
    }
    if (type2 < 0 || type2 > num_ssint_face_types - 1) {
        printf("Frog2 %d %d\n", type2, num_ssint_face_types - 1);
        return;
    }

    LJ_pair *p = &(params[LJI(type1, type2)]);
    *Emin = p->Emin;
    *Rmin = p->Rmin;
}

int LJ_matrix::get_num_types() {
    return num_ssint_face_types;
}

