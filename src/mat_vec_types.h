#ifndef MAT_VEC_TYPES_H_INCLUDED
#define MAT_VEC_TYPES_H_INCLUDED

/*
 * Defines what is meant by a scalar (essentially sets the precision of
 * the code between float or double).
 */
typedef double scalar;
//typedef long double scalar;

/*
 * A simple 3 dimensional vector (x, y, z)
 */
typedef struct
{
	scalar x, y, z;
} vector3;

/*
 * Defines a 12 vector (just for ease of use, clarity and compiler type checking)
 */
typedef scalar vector12[12];

/*
 * Defines a 12x12 matrix
 */
typedef scalar matrix12[12][12];

/*
 * Defines a 3x3 matrix
 */
typedef scalar matrix3[3][3];

/*
 * Defines a 4x4 matrix
 */
typedef scalar matrix4[4][4];

/*
 * A useful type for holding the upper triangular part of symmetric 4x4 matrices
 */
typedef struct
{
	scalar u00,	u01,	u02,	u03,
			u11,	u12,	u13,
				u22,	u23,
					u33;
} upper_triangular_matrix4;
#endif
