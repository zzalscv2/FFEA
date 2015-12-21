#include "mat_vec_fns_II.h"
#include <iostream>
using namespace std; 

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////
//
///  check whether two scalars have the same sign 
bool sameSign(scalar a, scalar b){
 return a < ffea_const::zero == b < ffea_const::zero;
}

/** Given 3 integers n0, n1, n2
 *    return the index missing in the list [0,1,2,3] 
 */
int getMissingNode(int n0, int n1, int n2) {

  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        return 3;
      } else {
        return 2;
      }
    } else {
      return 1;
    }
  } else {
    return 0;
  }

}



////////////////// SECTION 1 ///////////////////////
///  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////
/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
void arr3arr3Add(arr3 &vecA, arr3 &vecB, arr3 &res){

   for (int i=0; i<3; i++) {
     res[i] = vecA[i] + vecB[i];
    }

}


/** res = vecA - vecB */
/** res can be either vecA or vecB */
void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res){

   for (int i=0; i<3; i++) {
     res[i] = vecA[i] - vecB[i];
    }

}

/** w = u x v */
/**  (w != u) && (w != v) */
void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w)){

    w[0] = u[1]*v[2] - v[1]*u[2];
    w[1] = -u[0]*v[2] + v[0]*u[2];
    w[2] = u[0]*v[1] - v[0]*u[1];


}

/** return the dot product for arrays vecA and vecB */
scalar arr3arr3DotProduct(arr3 &vecA, arr3 &vecB) {

  scalar result = 0;
  for (int i=0; i<3; i++) {
     // cout << "vecA[" << i << "]: " << vecA[i] << " vecB[" << i << "]: " << vecB[i] << endl; 
     result += vecA[i] * vecB[i];
  }
  return result;
}


/** Normalise vector arr3 e */
void arr3Normalise(arr3 &e){

   scalar norm = 0.0;
   for (int i=0; i<3; i++) {
     norm += e[i]*e[i];
   }
   norm = sqrt(norm);
   for (int i=0; i<3; i++) {
     e[i] /= norm;
   }

}

/** get the normalised vector of arr3 e into arr3 n */
void arr3Normalise2(arr3 &e, arr3 &n){

   scalar norm = 0.0;
   for (int i=0; i<3; i++) {
     norm += e[i]*e[i];
   }
   norm = sqrt(norm);
   for (int i=0; i<3; i++) {
     n[i] = e[i]/norm;
   }

}


/** resize vector u, given scalar f */
void arr3Resize(scalar f, arr3 &u){

   for (int i=0; i<3; i++) {
      u[i] *= f;
    }

}

/** resize vector u into vector v, given scalar f */
void arr3Resize2(scalar f, arr3 &u, arr3 &v){

   for (int i=0; i<3; i++) {
      v[i] = f*u[i];
    }

}

/** Return the length of a vector v */
scalar mag(arr3 &v) {
   
   scalar s;
   for (int i=0; i<3; i++) {
      s += v[i] * v[i];
   } 
   return sqrt(s); 

}

void arr3Initialise(arr3 &v){ 
  
    for (int i=0; i<3; i++) {
      v[i] = ffea_const::zero; 
    }

} 
 



////////////////////////////////////////////////
////////////// END OF SECTION 1 ////////////////
////////////////////////////////////////////////



///////////////// SECTION 2 ////////////////////
////// Geometric functions for arr3 types //////
////////////////////////////////////////////////
/** t = unit(vecA - vecB) */
/** t can be either vecA or vecB */
void tangent(arr3 &vecA, arr3 &vecB, arr3 &t){

   scalar w=0;
   for (int i=0; i<3; i++) {
     t[i] = vecA[i] - vecB[i];
     w += t[i] * t[i];
    }
    w = sqrt(w);
    for (int i=0; i<3; i++) {
      t[i] /= w;
    }

}


/** w = unit(u x v) */
/**  (w != u) && (w != v) */ 
void getUnitNormal(arr3 &u, arr3 &v, arr3 &w){

   w[0] = u[1]*v[2] - v[1]*u[2];
   w[1] = -u[0]*v[2] + v[0]*u[2];
   w[2] = u[0]*v[1] - v[0]*u[1];
   scalar l = sqrt( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] );
   for (int i=0; i<3; i++) {
     w[i] /= l;
   }
}


/** calculate the normal vector n to the plane defined by the three points */
void getNormal(arr3 &v1, arr3 &v2, arr3 &v3, arr3 &n){

   arr3 pl1, pl2;
   arr3arr3Substract(v2, v1, pl1);
   arr3arr3Substract(v3, v1, pl2);
   getUnitNormal(pl1, pl2, n);

}

/* Given the face formed by tetA[0]:tetA[1]:tetA[2] 
 * get n, the normal to a face pointing inwards.
 */
void getNormalInwards(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &n){

   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   arr3 pl1, pl2, aux;
   arr3arr3Substract(tetA[n1], tetA[n0], pl1);
   arr3arr3Substract(tetA[n2], tetA[n0], pl2);
   // n is a unit vector normal to the face:
   arr3arr3VectorProduct(pl1, pl2, aux);
   arr3Normalise2(aux, n);
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   int n3 = getMissingNode(n0, n1, n2);
   arr3arr3Add(aux, tetA[n0], aux); 
   scalar d = - arr3arr3DotProduct(n, tetA[n0]);
   scalar t1 = arr3arr3DotProduct(tetA[n3], n) + d;
   scalar t2 = arr3arr3DotProduct(aux, n) + d;
   if (!sameSign(t1,t2)) arr3Resize(ffea_const::mOne, n);

}


/** check if points vec and test are at the same side
 *  of the plane formed by p1, p2 and p3 
 */
bool sameSidePlane(arr3 &vec, arr3 &test, arr3 &p1, arr3 &p2, arr3 &p3){

   arr3 pl1, pl2, n;
   arr3arr3Substract(p2, p1, pl1);
   arr3arr3Substract(p3, p1, pl2);
   arr3arr3VectorProduct(pl1, pl2, n);
   arr3Normalise(n); 
   scalar d = - arr3arr3DotProduct(n, p1);
   scalar t1 = arr3arr3DotProduct(vec, n) + d;
   scalar t2 = arr3arr3DotProduct(test, n) + d;
   // if ((t1 * t2) >= 0 ) return true;
   // cout << "t1: " << t1 << " t2: " << t2 << endl; 
   if (sameSign(t1,t2)) return true; 
   else return false;

}

/**
 */
bool samePlane(arr3 &p1, arr3 &p2, arr3 &p3, arr3 &p4){
 
   arr3 pl21, pl31, pl41; 
   arr3arr3Substract(p2, p1, pl21);
   arr3arr3Substract(p3, p1, pl31);
   arr3arr3Substract(p4, p1, pl41);
   scalar s = (pl21[1] * pl31[2] - pl21[2] * pl31[1])* pl41[0] + 
              (pl21[2] * pl31[0] - pl21[0] * pl31[2]) * pl41[1] + 
              (pl21[0] * pl31[1] - pl21[1] * pl31[0]) * pl41[2]; 
   if (fabs(s) > ffea_const::threeErr) return false; 
   else return true; 
   

}

/**  Given 4 co-planar points, check if ip and p1 lay on the same side 
 *     of the of the line formed by p2 and p3. 
 *   More specifically we check whether pl23 x pl21 and pl23 x pl2e 
 *     are parallel or antiparallel. 
 */
bool sameSideLine(arr3 &e, arr3 &p1, arr3 &p2, arr3 &p3) {

   // if (!samePlane(e, p1, p2, p3)) cout << "alarm" << endl;
   arr3 pl23, pl21, pl2e, v1, ve;
   arr3arr3Substract(p3, p2, pl23);
   arr3arr3Substract(p1, p2, pl21);
   arr3arr3Substract(e, p2, pl2e);
   arr3arr3VectorProduct(pl21, pl23, v1);
   arr3arr3VectorProduct(pl2e, pl23, ve);
   if (arr3arr3DotProduct(v1,ve) >= 0) return true;
   return false;

}

/**  check whether vector vec is in tetrahedron B. */
/**  more specifically, it will be there if 
 *     for each plane of the tetrahedron, 
 *     the point is on the same side as the remaining vertex */
bool nodeInTet(arr3 &vec, arr3 (tet)[4]){

   if (!sameSidePlane(vec, tet[0], tet[1], tet[2], tet[3])) return false;
   if (!sameSidePlane(vec, tet[1], tet[2], tet[3], tet[0])) return false;
   if (!sameSidePlane(vec, tet[2], tet[3], tet[0], tet[1])) return false;
   if (!sameSidePlane(vec, tet[3], tet[0], tet[1], tet[2])) return false;

   return true;

}

/** find the intersection point of the line that passes through the points e1 and e2, 
 *   and the plane defined by points p1, p2 and p3.
 *  \warning {this function should be called ONLY in the case 
 *            that intersection is known to occur.} 
 */
void linePlaneIntersectionPoint(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &p1, arr3 &p2, arr3 &p3) {

   // v is the vector of the line L(t) = e1 + v*t
   arr3 v;
   arr3arr3Substract(e2, e1, v);

   // now we need the unit vector that defines the plane, pn:
   arr3 pl1, pl2, pn;
   arr3arr3Substract(p2, p1, pl1);
   arr3arr3Substract(p3, p2, pl2);
   getUnitNormal(pl1, pl2, pn);

   // the plane is defined through: ax + by + cz + d = 0; 
   //   (a,b,c) = pn
   // so to find d we simply:
   scalar d = - arr3arr3DotProduct(pn, p1);

   // now find t and the point:
   scalar t = - (arr3arr3DotProduct(pn, e1) + d) / arr3arr3DotProduct(pn, v);
   arr3Resize(t, v);
   arr3arr3Add(e1, v, ip);

}

/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * more specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (tet[f2] - tet[f1]) and (tet[f3] - tet[f1]).
 *    - the intersection of a line is a point in the plane 
 */
bool intersectionPoint(arr3 &(ip), arr3 (&e1), arr3 (&e2), arr3 (&tet)[4], int f1, int f2, int f3){

  // check it e1 and e2 are on the same side of the plane:
  if ( sameSidePlane(e1, e2, tet[f1], tet[f2], tet[f3]) ) return false;

  // given that they are on different sides of the plane look for the intersection point.
  linePlaneIntersectionPoint(ip, e1, e2, tet[f1], tet[f2], tet[f3]);

  // and finally check whether this point ip belongs to the triangular face:
  if ( isPointInFace (ip, tet[f1], tet[f2], tet[f3]) ) return true;


  return false;

}


/** Check whether point ip is  
  *    inside of the three half-planes formed by the triangle's edges p1, p2, p3.
  */ 
bool isPointInFace(arr3 &ip, arr3 &p1, arr3 &p2, arr3 &p3) {

  if (! sameSideLine(ip, p1, p2, p3) ) return false;
  if (! sameSideLine(ip, p3, p1, p2) ) return false;
  if (! sameSideLine(ip, p2, p3, p1) ) return false;

  return true;
}

/** Return the center of coordinates for three points p1, p2, p3 in c */
void faceCentroid(arr3 &p1, arr3 &p2, arr3 &p3, arr3 &c){

   for (int i=0; i<3; i++){
     c[i] = ffea_const::oneOverThree * (p1[i] + p2[i] + p3[i]); 
   } 

}

///////////////// SECTION 3 ////////////////////
/// Transition functions from vector3 to arr3 // 
////////////////////////////////////////////////
void vec3Vec3SubsToArr3(vector3 &u, vector3 &v, arr3 (&w)){

    w[0] = u.x - v.x;
    w[1] = u.y - v.y;
    w[2] = u.z - v.z;

}

void vec3Arr3SubsToArr3(vector3 &u, arr3 (&v), arr3 (&w)){

    w[0] = u.x - v[0]; 
    w[1] = u.y - v[1]; 
    w[2] = u.z - v[2]; 

}

void vec3Arr3AddToArr3(vector3 &u, arr3 (&v), arr3 (&w)){

    w[0] = u.x + v[0]; 
    w[1] = u.y + v[1]; 
    w[2] = u.z + v[2]; 

}

void vec3ResizeToArr3(scalar f, vector3 &u, arr3 (&v)){

    v[0] = f*u.x;
    v[1] = f*u.y;
    v[2] = f*u.z;

}

scalar vec3Arr3DotProduct(vector3 &u, arr3 &v) {
    
     scalar s;
     s = u.x * v[0]; 
     s += u.y * v[1];
     s += u.z * v[2]; 
     return s; 

}


////////////// END OF SECTION 3 ////////////////



