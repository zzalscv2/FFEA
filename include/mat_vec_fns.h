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

#ifndef MAT_VEC_FNS_H_INCLUDED
#define MAT_VEC_FNS_H_INCLUDED

#include <math.h>
#include <stdio.h>

#include "mat_vec_types.h"
#include <Eigen/Geometry>

/**
 * Applies matrix A to vector v, storing result in v
 */
void mat12_apply(matrix12 A, vector12 v);

void vec3_mat3_mult(vector3 &v, matrix3 &A, vector3 &notV);
void vec3_mat3_mult_correct(vector3 &v, matrix3 &A, vector3 &notv);

void mat3_mult(matrix3 A, matrix3 B, matrix3 result);

void vec12_add(vector12 A, vector12 B);

void vec12_scale(vector12 A, scalar scale);

void mat3_mult_transpose(matrix3 A, matrix3 B, matrix3 result);

void mat3_mult_both_transposed(matrix3 A, matrix3 B, matrix3 result);

void mat3_scale(matrix3 A, scalar s);

scalar mat3_double_contraction_symmetric(matrix3 A);

scalar mat3_double_contraction(matrix3 A);

void get_rotation_matrix(vector3 &a, vector3 &b, matrix3 &R);

/** Inverts the given 3x3 matrix, storing the result in m_inv and the determinant in det_m */
void mat3_invert(matrix3 m, matrix3 m_inv, scalar *det_m);

void mat12_set_zero(matrix12 A);

void mat3_set_zero(matrix3 A);

void mat3_set_identity(matrix3 A);

void mat4_set_zero(matrix4 A);

void vector3_set_zero(vector3 &v);

void vec3_add_to_scaled(vector3 *v1, vector3 *v2, scalar a, int vec_size);

void vec3_scale_and_add(vector3 *v1, vector3 *v2, scalar a, int vec_size);

scalar vec3_normalise(vector3 &v1);

void vec12_set_zero(vector12 v);

/** Prints out the given 3x3 matrix */
void print_matrix3(matrix3 m);

/** Prints out the given 4x4 matrix */
void print_matrix4(matrix4 m);

/** Prints out the given 12x12 matrix */
void print_matrix12(matrix12 m);

/** Prints out the given 12-vector */
void print_vector12(vector12 v);

void print_vector3(vector3 &v);

#endif
