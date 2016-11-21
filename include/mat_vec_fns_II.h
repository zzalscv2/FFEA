#ifndef MAT_VEC_FNS_II_H_INCLUDED
#define MAT_VEC_FNS_II_H_INCLUDED

#include <cmath>
#include "mat_vec_types.h"

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////

/** check whether two scalars have the same sign */
template <class t_scalar> bool sameSign(t_scalar a, t_scalar b);

/** Given 3 integers n0, n1, n2
 *    return the index missing in the list [0,1,2,3] 
 */
int getMissingNode(int n0, int n1, int n2);

///////////////// SECTION 1 ////////////////////
////  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////

/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
template <class brr3> void arr3arr3Add(brr3 &vecA, brr3 &vecB, brr3 &res);

/** res = vecA - vecB */
/** res can be either vecA or vecB */
// void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res); 
template <class brr3> void arr3arr3Substract(brr3 &vecA, brr3 &vecB, brr3 &res);

/** w = u x v */
/** w must be different from u and v */
// void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w));
template <class brr3> void arr3arr3VectorProduct(brr3 (&u), brr3 (&v), brr3 (&w));

/** return the dot product for arrays vecA and vecB */
// scalar arr3arr3DotProduct(arr3 &vecA, arr3 &vecB);
template <class t_scalar, class brr3> t_scalar arr3arr3DotProduct(brr3 &vecA, brr3 &vecB);

/** Normalise vector arr3 e */
template <class t_scalar, class brr3> void arr3Normalise(brr3 &e);

/** Get the normalised vector of arr3 e into arr3 n */
// void arr3Normalise2(arr3 &e, arr3 &n);
template <class t_scalar, class brr3> void arr3Normalise2(brr3 &e, brr3 &n);

/** Given a scalar f, resize vector u*/
template <class t_scalar, class brr3> void arr3Resize(t_scalar f, brr3 &u);

/** Given a scalar f, resize vector u into vector v */
// void arr3Resize2(scalar f, arr3 &u, arr3 &v);
template <class t_scalar, class brr3> void arr3Resize2(t_scalar f, brr3 &u, brr3 &v);

template <class brr3> void arr3Store(brr3 &u, brr3 &v);

template <class t_scalar, class brr3> t_scalar arr3arr3Distance(brr3 &vecA, brr3 &vecB);

/** Return the length of a vector v */
template <class t_scalar> t_scalar mag(t_scalar (&v)[3]);

/** Initialise the input vector with (0, 0, 0) */
template <class brr3> void arr3Initialise(brr3 &v);


///////////////// SECTION 2 ////////////////////
///// Geometric  functions for arr3 types ////// 
////////////////////////////////////////////////
/** t = unit(vecA - vecB) */
/** t can be either vecA or vecB */
// void tangent(arr3 &vecA, arr3 &vecB, arr3 &t);
template <class t_scalar, class brr3> void tangent(brr3 &vecA, brr3 &vecB, brr3 &t);

/** w = unit(u x v) */
/**  (w != u) && (w != v) */
// template <class t_scalar, class brr3> void getUnitNormal(brr3 &u, brr3 &v, brr3 &w);
// template <class brr3> void getUnitNormal(brr3 &u, brr3 &v, brr3 &w);
template <class t_scalar, class brr3> void getUnitNormal(brr3 &u, brr3 &v, brr3 &w);

/** calculate the unit normal vector n to the plane defined by the three points */
template <class t_scalar, class brr3> void getNormal(brr3 &v1, brr3 &v2, brr3 &v3, brr3 &n);

/* Given the face formed by tetA[0]:tetA[1]:tetA[2] 
 * get n, the normal to a face pointing inwards.
 */
// void getNormalInwards(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &n);
template <class t_scalar, class brr3> void getNormalInwards(brr3 (&tetA)[4], int n0, int n1, int n2, brr3 &n);

/** check if points vec and test are at the same side
 *  of the plane formed by p1, p2 and p3 
 */
// bool sameSidePlane(arr3 &vec, arr3 &test, arr3 &p1, arr3 &p2, arr3 &p3);
template <class t_scalar, class brr3> bool sameSidePlane(brr3 &vec, brr3 &test, brr3 &p1, brr3 &p2, brr3 &p3);

/**  Given 4 co-planar points, check if ip and p1 lay on the same side 
 *     of the of the line formed by p2 and p3. 
 *   More specifically we check whether pl23 x pl21 and pl23 x pl2e 
 *     are parallel or antiparallel. 
 */
template <class t_scalar, class brr3> bool sameSideLine(brr3 &e, brr3 &p1, brr3 &p2, brr3 &p3);

/**  check whether vector vec is in tetrahedron B. */
/**  more specifically, it will be there if 
 *     for each plane of the tetrahedron, 
 *     the point is on the same side as the remaining vertex */
// bool nodeInTet(arr3 &vec, arr3 (tet)[4]);
template <class t_scalar, class brr3> bool nodeInTet(brr3 &vec, brr3 (tet)[4]);

/** Find the intersection point of the line that passes through the points e1 and e2, 
 *   and the plane defined by points p1, p2 and p3.
 *  \warning This function should be called ONLY in the case 
 *            that intersection is known to occur. 
 */
// void linePlaneIntersectionPoint(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &p1, arr3 &p2, arr3 &p3);
template <class t_scalar, class brr3> void linePlaneIntersectionPoint(brr3 &ip, brr3 &e1, brr3 &e2, brr3 &p1, brr3 &p2, brr3 &p3);

/** Check whether point ip is  
  *    inside of the three half-planes formed by the triangle's edges p1, p2, p3.
  */
// bool isPointInFace(arr3 &ip, arr3 &p1, arr3 &p2, arr3 &p3);
template <class t_scalar, class brr3> bool isPointInFace(brr3 &ip, brr3 &p1, brr3 &p2, brr3 &p3);

/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * More specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (tet[f2] - tet[f1]) and (tet[f3] - tet[f1]).
 *    - the intersection of a line is a point in the plane 
 */
// bool intersectionPoint(arr3 &(ip), arr3 (&e1), arr3 (&e2), arr3 (&tet)[4], int f1, int f2, int f3);
template <class t_scalar, class brr3> bool intersectionPoint(brr3 &(ip), brr3 (&e1), brr3 (&e2), brr3 (&tet)[4], int f1, int f2, int f3);

/** Return the center of coordinates for three points p1, p2, p3 in c */
// void faceCentroid(arr3 &p1, arr3 &p2, arr3 &p3, arr3 &c);
template <class brr3> void faceCentroid(brr3 &p1, brr3 &p2, brr3 &p3, brr3 &c);

/** Chech whether 4 points are on the same plane */
template <class t_scalar, class brr3> bool samePlane(brr3 &p1, brr3 &p2, brr3 &p3, brr3 &p4);

///////////////// SECTION 3 ////////////////////
//// Transition functions from vector3 to arr3 // 
////////////////////////////////////////////////
void vec3Vec3SubsToArr3(vector3 &u, vector3 &v, arr3 (&w));
void vec3Arr3SubsToArr3(vector3 &u, arr3 &v, arr3 &w);

void vec3Arr3AddToArr3(vector3 &u, arr3 (&v), arr3 (&w));
void vec3Vec3AddToArr3(vector3 &u, vector3 &v, arr3 (&w));

/** Given a scalar f, resize vec3 u into arr3 u*/
void vec3ResizeToArr3(scalar f, vector3 &u, arr3 &v);

scalar vec3Arr3DotProduct(vector3 &u, arr3 &v);

#endif 
