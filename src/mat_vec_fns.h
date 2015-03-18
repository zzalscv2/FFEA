#ifndef MAT_VEC_FNS_H_INCLUDED
#define MAT_VEC_FNS_H_INCLUDED

#include <math.h>

#include "mat_vec_types.h"

/*
 * Applies matrix A to vector v, storing result in v
 */
inline void mat12_apply(matrix12 A, vector12 v)
{
        int i, j;
        scalar temp_v[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

        for(i = 0; i < 12; i++)
                for(j = 0; j < 12; j++)
                        temp_v[i] += A[i][j]*v[j];
        for(i = 0; i < 12; i++) v[i] = temp_v[i];
}

/*
 *
 */
inline void mat3_mult(matrix3 A, matrix3 B, matrix3 result)
{
	int i, j, k;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
				result[i][j] += A[i][k] * B[k][j];
}

/*
 *
 */
inline void vec12_add(vector12 A, vector12 B) 
{
	for(int i = 0; i < 12; ++i) {
		B[i] = A[i] + B[i];
	}
}

/*
 *
 */
inline void mat3_mult_transpose(matrix3 A, matrix3 B, matrix3 result)
{
	int i, j, k;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
				result[i][j] += A[i][k] * B[j][k];
}

/*
 *
 */
inline void mat3_mult_both_transposed(matrix3 A, matrix3 B, matrix3 result)
{
	int i, j, k;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
				result[i][j] += A[k][i] * B[j][k];
}

/*
 *
 */
inline void mat3_scale(matrix3 A, scalar s)
{
	A[0][0] *= s; A[0][1] *= s; A[0][2] *= s;
	A[1][0] *= s; A[1][1] *= s; A[1][2] *= s;
	A[2][0] *= s; A[2][1] *= s; A[2][2] *= s;
}

/*
 *
 */
inline scalar mat3_double_contraction_symmetric(matrix3 A)
{
	return	A[0][0] * A[0][0] + A[1][1] * A[1][1] + A[2][2] * A[2][2]
		+ 2 * (A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2]);
}

/*
 *
 */
inline scalar mat3_double_contraction(matrix3 A)
{
	return	A[0][0] * A[0][0] + A[1][1] * A[1][1] + A[2][2] * A[2][2]
		+ A[0][1] * A[0][1] + A[1][0] * A[1][0] + A[0][2] * A[0][2]
		+ A[2][0] * A[2][0] + A[1][2] * A[1][2] + A[2][1] * A[2][1];
}

/*
 * Inverts the given 3x3 matrix, storing the result in m_inv and the determinant in det_m
 */
inline void mat3_invert(matrix3 m, matrix3 m_inv, scalar *det_m)
{
	scalar det;

	// Construct the inverse matrix
	m_inv[0][0] = m[2][2]*m[1][1] - m[2][1]*m[1][2];
	m_inv[0][1] = m[2][1]*m[0][2] - m[2][2]*m[0][1];
	m_inv[0][2] = m[1][2]*m[0][1] - m[1][1]*m[0][2];
	m_inv[1][0] = m[2][0]*m[1][2] - m[2][2]*m[1][0];
	m_inv[1][1] = m[2][2]*m[0][0] - m[2][0]*m[0][2];
	m_inv[1][2] = m[1][0]*m[0][2] - m[1][2]*m[0][0];
	m_inv[2][0] = m[2][1]*m[1][0] - m[2][0]*m[1][1];
	m_inv[2][1] = m[2][0]*m[0][1] - m[2][1]*m[0][0];
	m_inv[2][2] = m[1][1]*m[0][0] - m[1][0]*m[0][1];

	// calc determinant
	det = m[0][0] * m_inv[0][0] + m[1][0] * m_inv[0][1] + m[2][0] * m_inv[0][2];
	*det_m = det;

	// divide by determinant
	det = 1.0/det;
	m_inv[0][0]*=det; m_inv[0][1]*=det; m_inv[0][2]*=det;
	m_inv[1][0]*=det; m_inv[1][1]*=det; m_inv[1][2]*=det;
	m_inv[2][0]*=det; m_inv[2][1]*=det; m_inv[2][2]*=det;
}

/*
 *
 */
inline void mat12_set_zero(matrix12 A)
{
	int i, j;
	for(i = 0; i < 12; i++)
		for(j = 0; j < 12; j++)
			A[i][j] = 0;
}

/*
 *
 */
inline void mat3_set_zero(matrix3 A)
{
	int i, j;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			A[i][j] = 0;
}

/*
 *
 */
inline void mat4_set_zero(matrix4 A)
{
	int i, j;
	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			A[i][j] = 0;
}

/*
 *
 */
inline void vector3_set_zero(vector3 *v)
{
	v->x = 0; v->y = 0; v->z = 0;
}

/*
 *
 */ 
inline void vec3_add_to_scaled(vector3 *v1, vector3 *v2, scalar a, int vec_size)
{
	int i;
	#ifdef FFEA_PARALLEL_WITHIN_BLOB
		#pragma omp parallel for default(none) private(i) shared(v1, v2, a, vec_size)
	#endif
	for(i = 0; i < vec_size; i++) {
		v1[i].x += a * v2[i].x;
		v1[i].y += a * v2[i].y;
		v1[i].z += a * v2[i].z;
	}
}

/*
 *
 */ 
inline void vec3_scale_and_add(vector3 *v1, vector3 *v2, scalar a, int vec_size)
{
	int i;
	#ifdef FFEA_PARALLEL_WITHIN_BLOB
		#pragma omp parallel for default(none) private(i) shared(v1, v2, a, vec_size)
	#endif
	for(i = 0; i < vec_size; i++) {
		v1[i].x = a * v1[i].x + v2[i].x;
		v1[i].y = a * v1[i].y + v2[i].y;
		v1[i].z = a * v1[i].z + v2[i].z;
	}
}

/*
 *
 */
inline void vec12_set_zero(vector12 v)
{
	int i;
	for(i = 0; i < 12; i++)
		v[i] = 0;
}

/*
 * Prints out the given 3x3 matrix
 */
inline void print_matrix3(matrix3 m)
{
	int i;
	for(i = 0; i < 3; i++)
		printf("%e %e %e\n", m[i][0], m[i][1], m[i][2]);
}

/*
 * Prints out the given 4x4 matrix
 */
inline void print_matrix4(matrix4 m)
{
	int i;
	for(i = 0; i < 4; i++)
		printf("%e %e %e %e\n", m[i][0], m[i][1], m[i][2], m[i][3]);
}

/*
 * Prints out the given 12x12 matrix
 */
inline void print_matrix12(matrix12 m)
{
	int i, j;
	for(i = 0; i < 12; i++) {
		for(j = 0; j < 12; j++)
			printf("%e ", m[i][j]);
		printf("\n");
	}
}

/*
 * Prints out the given 12-vector
 */
inline void print_vector12(vector12 v)
{
	int i;
	for(i = 0; i < 12; i++)
		printf("%e\n", v[i]);
}

inline void print_vector3(vector3 *v)
{
	printf("%e %e %e\n", v->x, v->y, v->z);
}

inline double mag(vector3 *v)
{
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

inline vector3 normalise(vector3 *v) 
{
	double magnitude;
	vector3 norm;
	magnitude = mag(v);
	norm.x = v->x / magnitude;
	norm.y = v->y / magnitude;
	norm.z = v->z / magnitude;
	return norm;
}
#endif
