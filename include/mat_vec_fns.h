#ifndef MAT_VEC_FNS_H_INCLUDED
#define MAT_VEC_FNS_H_INCLUDED

#include <math.h>
#include <stdio.h>

#include "mat_vec_types.h"

/*
 * Applies matrix A to vector v, storing result in v
 */
void mat12_apply(matrix12 A, vector12 v);

void mat3_mult(matrix3 A, matrix3 B, matrix3 result);

void vec12_add(vector12 A, vector12 B);

void mat3_mult_transpose(matrix3 A, matrix3 B, matrix3 result);

void mat3_mult_both_transposed(matrix3 A, matrix3 B, matrix3 result);

void mat3_scale(matrix3 A, scalar s);

scalar mat3_double_contraction_symmetric(matrix3 A);

scalar mat3_double_contraction(matrix3 A);

/* Inverts the given 3x3 matrix, storing the result in m_inv and the determinant in det_m */
void mat3_invert(matrix3 m, matrix3 m_inv, scalar *det_m);

void mat12_set_zero(matrix12 A);

void mat3_set_zero(matrix3 A);

void mat4_set_zero(matrix4 A);

void vector3_set_zero(vector3 *v);

void vec3_add_to_scaled(vector3 *v1, vector3 *v2, scalar a, int vec_size);

void vec3_scale_and_add(vector3 *v1, vector3 *v2, scalar a, int vec_size);

void vec3_scale(vector3 *v, scalar scale);

void vec12_set_zero(vector12 v);

/* * Prints out the given 3x3 matrix */
void print_matrix3(matrix3 m);

/* * Prints out the given 4x4 matrix */
void print_matrix4(matrix4 m);

/* * Prints out the given 12x12 matrix */
void print_matrix12(matrix12 m);

/* * Prints out the given 12-vector */
void print_vector12(vector12 v);

void print_vector3(vector3 *v);

double mag(vector3 *v);

vector3 normalise(vector3 *v);

#endif
