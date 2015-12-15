#ifndef MAT_VEC_FNS_II_H_INCLUDED
#define MAT_VEC_FNS_II_H_INCLUDED

#include <cmath> 
#include "mat_vec_types.h"

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////
namespace ffea_const {
   const scalar mOne = -1.0;
   const scalar zero = 0.0;
   const scalar one = 1.0;
   const scalar two = 2.0;
   const scalar oneOverSix = 0.166666666666666667;
}

/** check whether two scalars have the same sign */
bool sameSign(scalar a, scalar b);

/** Given 3 integers n0, n1, n2
 *    return the index missing in the list [0,1,2,3] 
 */
int getMissingNode(int n0, int n1, int n2);


///////////////// SECTION 1 ////////////////////
///  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////

/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
void arr3arr3Add(arr3 &vecA, arr3 &vecB, arr3 &res);
 
/** res = vecA - vecB */
/** res can be either vecA or vecB */
void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res); 

/** w = u x v */
/** w must be different from u and v */
void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w));

/** return the dot product for arrays vecA and vecB */
scalar arr3arr3DotProduct(arr3 &vecA, arr3 &vecB);

/** Normalise vector arr3 e */
void arr3Normalise(arr3 &e);

/** get the normalised vector of arr3 e into arr3 n */
void arr3Normalise2(arr3 &e, arr3 &n);

/** Given a scalar f, resize vector u*/ 
void arr3Resize(scalar f, arr3 &u);


///////////////// SECTION 2 ////////////////////
///// Geometric  functions for arr3 types ////// 
////////////////////////////////////////////////
/** t = unit(vecA - vecB) */
/** t can be either vecA or vecB */
void tangent(arr3 &vecA, arr3 &vecB, arr3 &t);

/** w = unit(u x v) */
/**  (w != u) && (w != v) */
void getUnitNormal(arr3 &u, arr3 &v, arr3 &w);

/** calculate the unit normal vector n to the plane defined by the three points */
void getNormal(arr3 &v1, arr3 &v2, arr3 &v3, arr3 &n);

/* Given the face formed by tetA[0]:tetA[1]:tetA[2] 
 * get n, the normal to a face pointing inwards.
 */
void getNormalInwards(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &n);

/** check if points vec and test are at the same side
 *  of the plane formed by p1, p2 and p3 
 */
bool sameSidePlane(arr3 &vec, arr3 &test, arr3 &p1, arr3 &p2, arr3 &p3);

/**  Given 4 co-planar points, check if ip and p1 lay on the same side 
 *     of the of the line formed by p2 and p3. 
 *   More specifically we check whether pl23 x pl21 and pl23 x pl2e 
 *     are parallel or antiparallel. 
 */
bool sameSideLine(arr3 &e, arr3 &p1, arr3 &p2, arr3 &p3);


/**  check whether vector vec is in tetrahedron B. */
/**  more specifically, it will be there if 
 *     for each plane of the tetrahedron, 
 *     the point is on the same side as the remaining vertex */
bool nodeInTet(arr3 &vec, arr3 (tet)[4]);

/** find the intersection point of the line that passes through the points e1 and e2, 
 *   and the plane defined by points p1, p2 and p3.
 *  \warning {this function should be called ONLY in the case 
 *            that intersection is known to occur.} 
 */
void linePlaneIntersectionPoint(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &p1, arr3 &p2, arr3 &p3);

/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * more specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (tet[f2] - tet[f1]) and (tet[f3] - tet[f1]).
 *    - the intersection of a line is a point in the plane 
 */
bool intersectionPoint(arr3 &(ip), arr3 (&e1), arr3 (&e2), arr3 (&tet)[4], int f1, int f2, int f3);

/** Check whether point ip is  
  *    inside of the three half-planes formed by the triangle's edges p1, p2, p3.
  */ 
bool isPointInFace(arr3 &ip, arr3 &p1, arr3 &p2, arr3 &p3);



///////////////// SECTION 3 ////////////////////
/// Transition functions from vector3 to arr3 // 
////////////////////////////////////////////////
void vec3Vec3SubsToArr3(vector3 &u, vector3 &v, arr3 (&w));




#endif
