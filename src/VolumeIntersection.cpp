#include <iostream>
#include <cmath>
#include "VolumeIntersection.h"


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 * the code assumes that tetA has its vertices in the right order. 
 */
void getBAndN_Order(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &t, arr3 &b, arr3 &n){

  arr3  t_j; 
  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        getNormal(tetA[0], tetA[2], tetA[1], b);
      } else {
        getNormal(tetA[3], tetA[0], tetA[1], b);
      }
    } else {
        getNormal(tetA[3], tetA[2], tetA[0], b);
    }
  } else {
     getNormal(tetA[3], tetA[1], tetA[2], b);
  } 
  getUnitNormal(b, t, n); 
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  arr3 aux; 
  arr3arr3Add(tetA[n0], n, aux); 
  if (!sameSideLine(aux, tetA[n2], tetA[n1], tetA[n0])) arr3Resize(ffea_const::mOne, n); 
  // if we write the function we save a call to substract 
}


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 */
void getBAndN(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &t, arr3 &b, arr3 &n){
  
   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   //   but t == pl1, so: 
   arr3 pl2, aux; 
   arr3arr3Substract(tetA[n2], tetA[n0], pl2);
   // aux is a vector normal to the face:
   arr3arr3VectorProduct(t, pl2, aux);
   arr3Normalise2(aux, b);
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   int n3 = getMissingNode(n0, n1, n2); 
   scalar d = - arr3arr3DotProduct(b, tetA[n0]);
   scalar t1 = arr3arr3DotProduct(tetA[n3], b) + d;
   scalar t2 = arr3arr3DotProduct(aux, b) + d; // 1 + d; // arr3arr3DotProduct(b, b) + d   
   if (!sameSign(t1,t2)) arr3Resize(ffea_const::mOne, b);
 
  // Now we need N: 
  arr3arr3VectorProduct(b, t, aux); 
  arr3Normalise2(aux, n); 
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  if (!sameSideLine(aux, tetA[n2], tetA[n1], tetA[n0])) arr3Resize(ffea_const::mOne, n); 
  // if we were inlining the function we would save a call to substract 

}

/** Get the volume contribution of a node of type 1, i. e.,
 *    a node of the intersection polyhedron that already existed in one of the 
 *    intersectant thetrahedra. 
 */
scalar volumeForNode(arr3 (&tetA)[4], int node) {

  scalar volume = 0.0;

  // do the loop: 
  // for each edge:
  arr3 t_i, t_j, t_ij, n_j, b_j; 
  for (int i=0; i<4; i++) { 
    if (i == node) continue; 
    tangent(tetA[i], tetA[node], t_i); 
    scalar pt = arr3arr3DotProduct(tetA[node], t_i);
    // for each of the faces that this edge is adjacent to:
    for (int j=0; j<4; ++j) {
      if (j == node) continue;
      if (j == i) continue; 
      // getBAndN_Order(tetA, node, i, j, t_i, b_j, n_j); // both work well. 
      getBAndN(tetA, node, i, j, t_i, b_j, n_j);  // it works, too. 
      volume += pt * arr3arr3DotProduct(tetA[node], n_j) * arr3arr3DotProduct(tetA[node], b_j);
    } 
  } 

  return volume; 

}

/** Calculate the volume contribution of the intersection point ip,
 *    resulting of the intersection between edge tetA[e1]->tetA[e2]
 *    and the face given by the points tetB[f1]:tetB[f2]:tetB[f3]. 
 */ 
scalar volumeForIntPoint(arr3 &ip, arr3 (&tetA)[4], int e1, int e2, arr3 (&tetB)[4], int f1, int f2, int f3){

   scalar volume = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // arr3 B[3] will store the unit normals to these faces. 
   arr3 B[3]; 
   //
   // ip has three edges:
   //   tetA[e1]->tetA[e2],
   //   the intersection between face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   //   the intersection between face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   // arr3 T[3] will store the unit vectors of these edges.
   arr3 T[3]; 
   //
   // Finally, we will go through the loop edges x faces
   //
   /*
   cout << "intersection between edge tetA[" << e1 << "]:tetA[" << e2 << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
   cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;
   */ 
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards(tetB, f1, f2, f3, B[0]); 
   //  and the easy tangent:  
   tangent(tetA[e2], tetA[e1], T[0]); 
 
   // Secondly, normals and tangents for the tetA faces: 
   int cnt = 0; 
   int F[2];
   for (int i=0; i<4; i++){ 
     if ((i != e1) && (i != e2)) {
       // store "i" so we'll know that B[cnt] will correspond to the face triad e1:e2:i
       F[cnt] = i;
       cnt++; 
       getNormalInwards(tetA, e1, e2, i, B[cnt]);
       getUnitNormal(B[0],B[cnt],T[cnt]);
       // but obviously T[cnt] may not be in the right sense: 
       //   we know the intersection point, 
       //           and that it comes from the edge tetA[e1]:tetA[e2]
       //   Hence, for
       //   the intersection between face ( tetA[e1]->tetA[i] ) : (tetA[e1]->tetA[e2]) 
       //       and face tetB[f1]:tetB[f2]:tetB[f3];
       //   T[1] has to point towards... 
       int e4 = getMissingNode(i, e1, e2);
       arr3 vaux;
       arr3arr3Add(ip, T[cnt], vaux);
       if (!sameSidePlane( vaux , tetA[i], tetA[e4], tetA[e1], tetA[e2])) arr3Resize(ffea_const::mOne, T[cnt]);
     }
   }

  
   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., T[0] with faces B[1] and B[2]
   arr3 n; 
   scalar ipT[3], ipB[3];
   for (int i=0; i<3; i++) {
      ipT[i] = arr3arr3DotProduct(ip,T[i]); 
      ipB[i] = arr3arr3DotProduct(ip,B[i]); 
   } 
   for (int i=0; i<2; i++) { 
      // for B[1] we need N
      getUnitNormal(B[i+1], T[0], n); 
      // and we'll check whether n lays on the same side as tetA[F[i]]
      //    of the line given by T[0] = unit(tetA[e2] - tetA[e1])
      arr3 vaux;
      arr3arr3Add(ip, n, vaux);
      if (!sameSideLine(vaux, tetA[F[i]], tetA[e2], tetA[e1])) arr3Resize(ffea_const::mOne, n); 
      volume += ipT[0] * arr3arr3DotProduct(ip, n) * ipB[i+1]; //  arr3arr3DotProduct(ip, B[i+1]);
   } 
  
   // Now we need the contribution of the edges T[1] and T[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   arr3 V[3];
   arr3arr3Add(ip, T[1], V[0]);
   arr3arr3Add(ip, T[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First T[i+1] with B[0];
      // get the normal:
      getUnitNormal(B[0], T[i+1], n); 
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + T[i]*t  
      arr3 v;
      arr3arr3Add(ip, n, v);
      if (!sameSideLine(v, V[i+1], V[i], ip)) arr3Resize(ffea_const::mOne, n); 
      // volume += arr3arr3DotProduct(ip, T[i+1]) * arr3arr3DotProduct(ip, n) * arr3arr3DotProduct(ip, B[0]);
      volume += ipT[i+1] * arr3arr3DotProduct(ip, n) * ipB[0];   

      // and now T[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal(B[i+1],T[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (arr3arr3DotProduct(n,B[0]) < 0) arr3Resize(ffea_const::mOne, n); 
      volume += ipT[i+1] * arr3arr3DotProduct(ip, n) * ipB[i+1];
   } 
      
   return volume;
}


/** Return the volume intersection between two tetrahedra */ 
scalar volumeIntersection(arr3 (&tetA)[4], arr3 (&tetB)[4]){

  scalar vol = 0.0; 
  scalar aux;
  
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (nodeInTet(tetA[i], tetB)) { 
      vol += volumeForNode(tetA, i);
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (nodeInTet(tetB[i], tetA)) { 
      vol += volumeForNode(tetB, i);
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  arr3 ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        vol += volumeForIntPoint(ip, tetA, i, j, tetB, 0, 1, 2); 
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        vol += volumeForIntPoint(ip, tetA, i, j, tetB, 1, 2, 3); 
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        vol += volumeForIntPoint(ip, tetA, i, j, tetB, 2, 3, 0); 
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        vol += volumeForIntPoint(ip, tetA, i, j, tetB, 3, 0, 1); 
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  return vol; 
 

} 
