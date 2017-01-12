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

#include <iostream>
#include <cmath>
#include <iomanip>
#include "VolumeIntersection.h"

using namespace std; 

template <class t_scalar, class brr3> bool exists(brr3 &p, int ips, brr3 (&W)[56]){

  for (int i=0; i<ips; i++){ 
    if (arr3arr3Distance<t_scalar,brr3>(p, W[i]) < ffea_const::threeErr) {
      // cout << "discarded" << endl;
      return true;
    }
  } 
  return false;

}

/**  max volume for a number of points computed as the size of the cube 
 *    with the same side size as the double of the maximum distance found 
 *    between the given points.
 */
template <class t_scalar,class brr3> t_scalar maxVolume(int ips, brr3 (&W)[56]){

  brr3 c; 
  arr3Initialise(c); 
  for (int i=0; i<ips; i++){
    arr3arr3Add<t_scalar,brr3>(c, W[i], c); 
  } 
  t_scalar s = 1./ips; 
  arr3Resize<t_scalar,brr3>(s, c);
  t_scalar d0 = ffea_const::zero;
  t_scalar di; 
  for (int i=0; i<ips; i++){
    di = arr3arr3Distance<t_scalar,brr3>(c, W[i]); 
    if (di > d0) d0 = di;
  } 
  return ffea_const::eight*di*di*di; 
  // return ffea_const::volFactor*di*di*di; 

}

/**  max volume and area for a number of points computed as the size of the cube 
 *    with the same side size as the double of the maximum distance found 
 *    between the given points.
 */
template <class t_scalar,class brr3> void maxVolumeAndArea(int ips, brr3 (&W)[56], t_scalar &volume, t_scalar &area){

  brr3 c; 
  arr3Initialise(c); 
  for (int i=0; i<ips; i++){
    arr3arr3Add<t_scalar,brr3>(c, W[i], c); 
  } 
  t_scalar s = 1./ips; 
  arr3Resize<t_scalar,brr3>(s, c);
  t_scalar d0 = ffea_const::zero;
  t_scalar di; 
  for (int i=0; i<ips; i++){
    di = arr3arr3Distance<t_scalar,brr3>(c, W[i]); 
    if (di > d0) d0 = di;
  } 
  area = ffea_const::twentyfour*di*di;
  volume = ffea_const::eight*di*di*di; 
  // return ffea_const::volFactor*di*di*di; 

}


/**  CM for the first ips points stored in brr3 W[56] */ 
template <class t_scalar,class brr3> void findCM(int ips, brr3 (&W)[56], brr3 &cm){

  arr3Initialise(cm); 
  for (int i=0; i<ips; i++){
    arr3arr3Add<t_scalar,brr3>(cm, W[i], cm); 
  } 
  t_scalar s = 1./ips; 
  arr3Resize<t_scalar,brr3>(s, cm);

}


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 * the code assumes that tetA has its vertices in the right order. 
 */
template <class t_scalar,class brr3> void getBAndN_Order(brr3 (&tetA)[4], int n0, int n1, int n2, brr3 &t, brr3 &b, brr3 &n){

  // brr3  t_j; 
  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        getNormal<t_scalar,brr3>(tetA[0], tetA[2], tetA[1], b);
      } else {
        getNormal<t_scalar,brr3>(tetA[3], tetA[0], tetA[1], b);
      }
    } else {
        getNormal<t_scalar,brr3>(tetA[3], tetA[2], tetA[0], b);
    }
  } else {
     getNormal<t_scalar,brr3>(tetA[3], tetA[1], tetA[2], b);
  } 
  getUnitNormal<t_scalar,brr3>(b, t, n); 
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  brr3 aux; 
  arr3arr3Add<t_scalar,brr3>(tetA[n0], n, aux); 
  if (!sameSideLine<t_scalar,brr3>(aux, tetA[n2], tetA[n1], tetA[n0])) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
  // if we write the function we save a call to substract 
}


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 */
// void getBAndN(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &t, arr3 &b, arr3 &n){
template <class t_scalar, class brr3> void getBAndN(brr3 (&tetA)[4], int n0, int n1, int n2, brr3 &t, brr3 &b, brr3 &n){
  
   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   //   but t == pl1, so: 
   brr3 pl2, aux; 
   arr3arr3Substract<t_scalar,brr3>(tetA[n2], tetA[n0], pl2);
   // aux is a vector normal to the face:
   arr3arr3VectorProduct<t_scalar,brr3>(t, pl2, aux);
   arr3Normalise2<t_scalar,brr3>(aux, b);
   arr3arr3Add<t_scalar,brr3>(aux, tetA[n0], aux); 
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   int n3 = getMissingNode(n0, n1, n2); 
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(b, tetA[n0]);
   t_scalar t1 = arr3arr3DotProduct<t_scalar,brr3>(tetA[n3], b) + d;
   t_scalar t2 = arr3arr3DotProduct<t_scalar,brr3>(aux, b) + d; // 1 + d; // arr3arr3DotProduct(b, b) + d   
   if (!sameSign(t1,t2)) arr3Resize<t_scalar,brr3>(ffea_const::mOne, b);
 
  // Now we need N: 
  getUnitNormal<t_scalar,brr3>(b, t, n);
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  arr3arr3Add<t_scalar,brr3>(tetA[n0], n, aux);
  if (!sameSideLine<t_scalar,brr3>(aux, tetA[n2], tetA[n1], tetA[n0])) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n);

}

/** Get the volume contribution of a node of type 1, i. e.,
 *    a node of the intersection polyhedron that already existed in one of the 
 *    intersectant thetrahedra. 
 */
template <class t_scalar, class brr3> t_scalar volumeForNode(brr3 (&tetA)[4], int node) {

  t_scalar volume = 0.0;

  // do the loop: 
  // for each edge:
  brr3 t_i, n_j, b_j; // t_j, t_ij
  for (int i=0; i<4; i++) { 
    if (i == node) continue; 
    tangent<t_scalar,brr3>(tetA[i], tetA[node], t_i); 
    t_scalar pt = arr3arr3DotProduct<t_scalar,brr3>(tetA[node], t_i);
    // for each of the faces that this edge is adjacent to:
    for (int j=0; j<4; ++j) {
      if (j == node) continue;
      if (j == i) continue; 
      // getBAndN_Order<t_scalar,brr3>(tetA, node, i, j, t_i, b_j, n_j); // both work well. 
      getBAndN<t_scalar,brr3>(tetA, node, i, j, t_i, b_j, n_j);  // it works, too. 
      /////////////////
      /* // CHECK ! // 
      brr3 aux; 
      cout << "node: " << node << " i: " << i << " j: " << j << endl; 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,t_i,aux); 
      arr3arr3Add<t_scalar,brr3>(tetA[node], aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,b_j,aux); 
      arr3arr3Add<t_scalar,brr3>(tetA[node], aux, aux); 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,n_j,aux); 
      arr3arr3Add<t_scalar,brr3>(tetA[node], aux, aux); 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! // 
      /////////////////
      volume += pt * arr3arr3DotProduct<t_scalar,brr3>(tetA[node], n_j) * arr3arr3DotProduct<t_scalar,brr3>(tetA[node], b_j);
    } 
  } 

  // cout << " volume: " << volume << endl; 
  return volume; 

}

/** Get the volume and area contribution of a node of type 1, i. e.,
 *    a node of the intersection polyhedron that already existed in one of the 
 *    intersectant thetrahedra. 
 */
template <class t_scalar, class brr3> void volumeAndAreaForNode(brr3 (&tetA)[4], int node, t_scalar &volume, t_scalar &area) {

  volume = 0.0;
  area = 0.0; 

  // do the loop: 
  // for each edge:
  brr3 t_i, n_j, b_j; // t_j, t_ij
  for (int i=0; i<4; i++) { 
    if (i == node) continue; 
    tangent<t_scalar,brr3>(tetA[i], tetA[node], t_i); 
    t_scalar pt = arr3arr3DotProduct<t_scalar,brr3>(tetA[node], t_i);
    // for each of the faces that this edge is adjacent to:
    for (int j=0; j<4; ++j) {
      if (j == node) continue;
      if (j == i) continue; 
      // getBAndN_Order<t_scalar,brr3>(tetA, node, i, j, t_i, b_j, n_j); // both work well. 
      getBAndN<t_scalar,brr3>(tetA, node, i, j, t_i, b_j, n_j);  // it works, too. 
      t_scalar pt_pn = pt * arr3arr3DotProduct<t_scalar,brr3>(tetA[node], n_j);
      area += pt_pn; 
      volume += pt_pn * arr3arr3DotProduct<t_scalar,brr3>(tetA[node], b_j);
    } 
  } 
}


/** Calculate the volume contribution of the intersection point ip,
 *    resulting of the intersection between edge tetA[e1]->tetA[e2]
 *    and the face given by the points tetB[f1]:tetB[f2]:tetB[f3]. 
 */ 
template <class t_scalar, class brr3> t_scalar volumeForIntPoint(brr3 &ip, brr3 (&tetA)[4], int e1, int e2, brr3 (&tetB)[4], int f1, int f2, int f3){

   t_scalar volume = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // arr3 B[3] will store the unit normals to these faces. 
   brr3 B[3]; 
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
   brr3 T[3]; 
   //
   // Finally, we will go through the loop edges x faces. 
   //
   // Let's start:
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards<t_scalar,brr3>(tetB, f1, f2, f3, B[0]); 
   //  and the easy tangent:  
   tangent<t_scalar,brr3>(tetA[e2], tetA[e1], T[0]); 
   // but check whether we are entering or escaping: 
   int f4 = getMissingNode(f1,f2,f3); 
   //   so if escaping:
   if (sameSidePlane<t_scalar,brr3>(tetA[e1], tetB[f4], tetB[f1], tetB[f2], tetB[f3])) 
     arr3Resize<t_scalar,brr3>(ffea_const::mOne, T[0]); 
 
   // Secondly, normals and tangents for the tetA faces: 
   int cnt = 0; 
   int F[2];
   for (int i=0; i<4; i++){ 
     if ((i != e1) && (i != e2)) {
       // store "i" so we'll know that B[cnt] will correspond to the face triad e1:e2:i
       F[cnt] = i;
       cnt++; 
       getNormalInwards<t_scalar,brr3>(tetA, e1, e2, i, B[cnt]);
       getUnitNormal<t_scalar,brr3>(B[0],B[cnt],T[cnt]);
       // but obviously T[cnt] may not be in the right sense: 
       //   we know the intersection point, 
       //           and that it comes from the edge tetA[e1]:tetA[e2]
       //   Hence, for
       //   the intersection between face ( tetA[e1]->tetA[i] ) : (tetA[e1]->tetA[e2]) 
       //       and face tetB[f1]:tetB[f2]:tetB[f3];
       //   T[1] has to point towards... 
       int e4 = getMissingNode(i, e1, e2);
       brr3 vaux;
       arr3arr3Add<t_scalar,brr3>(ip, T[cnt], vaux);
       if (!sameSidePlane<t_scalar,brr3>( vaux , tetA[i], tetA[e4], tetA[e1], tetA[e2])) 
                 arr3Resize<t_scalar,brr3>(ffea_const::mOne, T[cnt]);
     }
   }

  /////////////////////
  ////// CHECK ////////
  /*
  brr3 aux; 
  brr3 C[3]; 
  faceCentroid(tetB[f1], tetB[f2], tetB[f3], C[0]);
  faceCentroid(tetA[e1], tetA[e2], tetA[F[0]], C[1]);
  faceCentroid(tetA[e1], tetA[e2], tetA[F[1]], C[2]);
  arr3arr3Add<t_scalar,brr3>(C[0], B[0], aux); 
  cout << "C: " << C[0][0] << ", " << C[0][1] << ", " << C[0][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 

  arr3arr3Add<t_scalar,brr3>(C[1], B[1], aux); 
  cout << "C: " << C[1][0] << ", " << C[1][1] << ", " << C[1][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 
  arr3arr3Add<t_scalar,brr3>(C[2], B[2], aux); 
  cout << "C: " << C[2][0] << ", " << C[2][1] << ", " << C[2][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 
  */ 
  ///////////////////////////
  /*
  for (int i=0; i<3; i++) { 
    arr3arr3Add<t_scalar,brr3>(ip, T[i], aux); 
    cout << "P: " << ip[0] << " " << ip[1] << " " << ip[2] << endl; 
    cout << "T: " << aux[0] << " " << aux[1] << " " << aux[2] << endl; 
    cout << endl << endl; 
  } 
  */ 
  ////// CHECK ////////
  /////////////////////
  
   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., T[0] with faces B[1] and B[2]
   brr3 n; 
   t_scalar ipT[3], ipB[3];
   for (int i=0; i<3; i++) {
      ipT[i] = arr3arr3DotProduct<t_scalar,brr3>(ip,T[i]); 
      ipB[i] = arr3arr3DotProduct<t_scalar,brr3>(ip,B[i]); 
   } 
   for (int i=0; i<2; i++) { 
      // for B[1] we need N
      getUnitNormal<t_scalar,brr3>(B[i+1], T[0], n); 
      // and we'll check whether n lays on the same side as tetA[F[i]]
      //    of the line given by T[0] = unit(tetA[e2] - tetA[e1])
      brr3 vaux;
      arr3arr3Add<t_scalar,brr3>(ip, n, vaux);
      if (!sameSideLine<t_scalar,brr3>(vaux, tetA[F[i]], tetA[e2], tetA[e1])) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      /////////////////
      /*// CHECK ! // 
      brr3 aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,T[0],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,B[i+1],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,n,aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      // CHECK ! // */
      /////////////////
      volume += ipT[0] * arr3arr3DotProduct<t_scalar,brr3>(ip, n) * ipB[i+1]; //  arr3arr3DotProduct(ip, B[i+1]);
   } 
  
   // Now we need the contribution of the edges T[1] and T[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   brr3 V[3];
   arr3arr3Add<t_scalar,brr3>(ip, T[1], V[0]);
   arr3arr3Add<t_scalar,brr3>(ip, T[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First T[i+1] with B[0];
      // get the normal:
      getUnitNormal<t_scalar,brr3>(B[0], T[i+1], n); 
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + T[i]*t  
      brr3 v;
      arr3arr3Add<t_scalar,brr3>(ip, n, v);
      if (!sameSideLine<t_scalar,brr3>(v, V[i+1], V[i], ip)) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      volume += ipT[i+1] * arr3arr3DotProduct<t_scalar,brr3>(ip, n) * ipB[0];   

      /////////////////
      /*// CHECK ! // 
      brr3 aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,T[i+1],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,B[0],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,n,aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! //
      /////////////////

      // and now T[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal<t_scalar,brr3>(B[i+1],T[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (arr3arr3DotProduct<t_scalar,brr3>(n,B[0]) < 0) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      volume += ipT[i+1] * arr3arr3DotProduct<t_scalar,brr3>(ip, n) * ipB[i+1];
      /////////////////
      /* // CHECK ! // 
      brr3 aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,T[i+1],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,B[i+1],aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      arr3Resize2<t_scalar,brr3>(ffea_const::ten,n,aux); 
      arr3arr3Add<t_scalar,brr3>(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! //
      /////////////////

   } 
      
   // cout << " volume: " << volume << endl; 
   return volume;
}


/** Calculate the volume and area contribution of the intersection point ip,
 *    resulting of the intersection between edge tetA[e1]->tetA[e2]
 *    and the face given by the points tetB[f1]:tetB[f2]:tetB[f3]. 
 */ 
template <class t_scalar, class brr3> void volumeAndAreaForIntPoint(brr3 &ip, brr3 (&tetA)[4], int e1, int e2, brr3 (&tetB)[4], int f1, int f2, int f3, t_scalar &volume, t_scalar &area){

   volume = 0; 
   area = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // arr3 B[3] will store the unit normals to these faces. 
   brr3 B[3]; 
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
   brr3 T[3]; 
   //
   // Finally, we will go through the loop edges x faces. 
   //
   // Let's start:
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards<t_scalar,brr3>(tetB, f1, f2, f3, B[0]); 
   //  and the easy tangent:  
   tangent<t_scalar,brr3>(tetA[e2], tetA[e1], T[0]); 
   // but check whether we are entering or escaping: 
   int f4 = getMissingNode(f1,f2,f3); 
   //   so if escaping:
   if (sameSidePlane<t_scalar,brr3>(tetA[e1], tetB[f4], tetB[f1], tetB[f2], tetB[f3])) 
     arr3Resize<t_scalar,brr3>(ffea_const::mOne, T[0]); 
 
   // Secondly, normals and tangents for the tetA faces: 
   int cnt = 0; 
   int F[2];
   for (int i=0; i<4; i++){ 
     if ((i != e1) && (i != e2)) {
       // store "i" so we'll know that B[cnt] will correspond to the face triad e1:e2:i
       F[cnt] = i;
       cnt++; 
       getNormalInwards<t_scalar,brr3>(tetA, e1, e2, i, B[cnt]);
       getUnitNormal<t_scalar,brr3>(B[0],B[cnt],T[cnt]);
       // but obviously T[cnt] may not be in the right sense: 
       //   we know the intersection point, 
       //           and that it comes from the edge tetA[e1]:tetA[e2]
       //   Hence, for
       //   the intersection between face ( tetA[e1]->tetA[i] ) : (tetA[e1]->tetA[e2]) 
       //       and face tetB[f1]:tetB[f2]:tetB[f3];
       //   T[1] has to point towards... 
       int e4 = getMissingNode(i, e1, e2);
       brr3 vaux;
       arr3arr3Add<t_scalar,brr3>(ip, T[cnt], vaux);
       if (!sameSidePlane<t_scalar,brr3>( vaux , tetA[i], tetA[e4], tetA[e1], tetA[e2])) 
                 arr3Resize<t_scalar,brr3>(ffea_const::mOne, T[cnt]);
     }
   }

   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., T[0] with faces B[1] and B[2]
   brr3 n; 
   t_scalar ipT[3], ipB[3];
   for (int i=0; i<3; i++) {
      ipT[i] = arr3arr3DotProduct<t_scalar,brr3>(ip,T[i]); 
      ipB[i] = arr3arr3DotProduct<t_scalar,brr3>(ip,B[i]); 
   } 
   for (int i=0; i<2; i++) { 
      // for B[1] we need N
      getUnitNormal<t_scalar,brr3>(B[i+1], T[0], n); 
      // and we'll check whether n lays on the same side as tetA[F[i]]
      //    of the line given by T[0] = unit(tetA[e2] - tetA[e1])
      brr3 vaux;
      arr3arr3Add<t_scalar,brr3>(ip, n, vaux);
      if (!sameSideLine<t_scalar,brr3>(vaux, tetA[F[i]], tetA[e2], tetA[e1])) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      t_scalar ipT_ipN = ipT[0] * arr3arr3DotProduct<t_scalar,brr3>(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[i+1];
   } 
  
   // Now we need the contribution of the edges T[1] and T[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   brr3 V[3];
   arr3arr3Add<t_scalar,brr3>(ip, T[1], V[0]);
   arr3arr3Add<t_scalar,brr3>(ip, T[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First T[i+1] with B[0];
      // get the normal:
      getUnitNormal<t_scalar,brr3>(B[0], T[i+1], n); 
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + T[i]*t  
      brr3 v;
      arr3arr3Add<t_scalar,brr3>(ip, n, v);
      if (!sameSideLine<t_scalar,brr3>(v, V[i+1], V[i], ip)) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      t_scalar ipT_ipN = ipT[i+1] * arr3arr3DotProduct<t_scalar,brr3>(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[0];   

      // and now T[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal<t_scalar,brr3>(B[i+1],T[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (arr3arr3DotProduct<t_scalar,brr3>(n,B[0]) < 0) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n); 
      ipT_ipN = ipT[i+1] * arr3arr3DotProduct<t_scalar,brr3>(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[i+1];
   } 
}


/** Return the volume intersection between two tetrahedra */ 
// scalar volumeIntersection(arr3 (&tetA)[4], arr3 (&tetB)[4]){
template <class t_scalar, class brr3> t_scalar volumeIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4]){ 


  t_scalar vol = 0.0; 
  brr3 W[56]; 

  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists<t_scalar,brr3>(tetA[i], ips, W) && nodeInTet<t_scalar,brr3>(tetA[i], tetB)) { 
      // cout << "tetA-node[" << i << "] found int tetB " << endl; 
      arr3Store<t_scalar,brr3>(tetA[i], W[ips]);
      ips += 1; 
      vol += volumeForNode<t_scalar,brr3>(tetA, i);
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists<t_scalar,brr3>(tetB[i], ips, W) && nodeInTet<t_scalar,brr3>(tetB[i], tetA)) { 
      // cout << "tetB-node[" << i << "] found int tetA " << endl; 
      arr3Store<t_scalar,brr3>(tetB[i], W[ips]);
      ips += 1;
      vol += volumeForNode<t_scalar,brr3>(tetB, i);
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  brr3 ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /* cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 0, 1, 2); 
        } 
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 1, 2, 3); 
        } 
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 2, 3, 0); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 3, 0, 1); 
        }
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 0, 1, 2); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 1, 2, 3); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 2, 3, 0); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 3, 0, 1); 
        }
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) 
    return 0.; 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  t_scalar w = maxVolume<t_scalar,brr3>(ips, W); 
  if ((fabs(vol) > w) || (vol < 0)) { 
    vol = w;
  }

  return vol; 
 

} 

/** Return the volume intersection between two tetrahedra */ 
template <class t_scalar, class brr3> t_scalar volumeIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4], brr3 &cm){


  t_scalar vol = 0.0; 
  brr3 W[56]; 

  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists<t_scalar,brr3>(tetA[i], ips, W) && nodeInTet<t_scalar,brr3>(tetA[i], tetB)) { 
      // cout << "tetA-node[" << i << "] found int tetB " << endl; 
      arr3Store<t_scalar,brr3>(tetA[i], W[ips]);
      ips += 1; 
      vol += volumeForNode<t_scalar,brr3>(tetA, i);
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists<t_scalar,brr3>(tetB[i], ips, W) && nodeInTet<t_scalar,brr3>(tetB[i], tetA)) { 
      // cout << "tetB-node[" << i << "] found int tetA " << endl; 
      arr3Store<t_scalar,brr3>(tetB[i], W[ips]);
      ips += 1;
      vol += volumeForNode<t_scalar,brr3>(tetB, i);
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  brr3 ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /* cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 0, 1, 2); 
        } 
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 1, 2, 3); 
        } 
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 2, 3, 0); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 3, 0, 1); 
        }
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 0, 1, 2); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 1, 2, 3); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 2, 3, 0); 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 3, 0, 1); 
        }
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) 
    return 0.; 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  t_scalar w = maxVolume<t_scalar,brr3>(ips, W); 
  if ((fabs(vol) > w) || (vol < 0)) { 
    vol = w;
  }

  findCM<t_scalar,brr3>(ips, W, cm);

  return vol; 
 

} 

/** Return the volume intersection between two tetrahedra */ 
template <class t_scalar, class brr3> t_scalar volumeIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4], brr3 &cm, brr3 &dr){


  t_scalar vol = 0.0; 
  brr3 W[56]; 
  brr3 cmA, cmB;
  brr3 tetAin[4], tetBin[4];
  int countAin[4], countBin[4];
  for (int i=0; i<4; i++){
    countAin[i] = 0;
    countBin[i] = 0;
  }


  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists<t_scalar,brr3>(tetA[i], ips, W) && nodeInTet<t_scalar,brr3>(tetA[i], tetB)) { 
      // cout << "tetA-node[" << i << "] found int tetB " << endl; 
      arr3Store<t_scalar,brr3>(tetA[i], W[ips]);
      ips += 1; 
      vol += volumeForNode<t_scalar,brr3>(tetA, i);
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists<t_scalar,brr3>(tetB[i], ips, W) && nodeInTet<t_scalar,brr3>(tetB[i], tetA)) { 
      // cout << "tetB-node[" << i << "] found int tetA " << endl; 
      arr3Store<t_scalar,brr3>(tetB[i], W[ips]);
      ips += 1;
      vol += volumeForNode<t_scalar,brr3>(tetB, i);
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  brr3 ip; 
  
  getNormalInwards<t_scalar,brr3>(tetB, 0, 1, 2, tetBin[0]); 
  getNormalInwards<t_scalar,brr3>(tetB, 1, 2, 3, tetBin[1]); 
  getNormalInwards<t_scalar,brr3>(tetB, 2, 3, 0, tetBin[2]); 
  getNormalInwards<t_scalar,brr3>(tetB, 3, 0, 1, tetBin[3]); 
  getNormalInwards<t_scalar,brr3>(tetA, 0, 1, 2, tetAin[0]); 
  getNormalInwards<t_scalar,brr3>(tetA, 1, 2, 3, tetAin[1]); 
  getNormalInwards<t_scalar,brr3>(tetA, 2, 3, 0, tetAin[2]); 
  getNormalInwards<t_scalar,brr3>(tetA, 3, 0, 1, tetAin[3]); 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /* cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 0, 1, 2);
        } 
        countBin[0] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 1, 2, 3); 
        } 
        countBin[1] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 2, 3, 0); 
        }
        countBin[2] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 3, 0, 1); 
        }
        countBin[3] += 1;
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 0, 1, 2); 
        }
        countAin[0] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 1, 2, 3); 
        }
        countAin[1] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 2, 3, 0); 
        }
        countAin[2] += 1;
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 3, 0, 1); 
        }
        countAin[3] += 1;
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) 
    return 0.; 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  t_scalar w = maxVolume<t_scalar,brr3>(ips, W); 
  if ((fabs(vol) > w) || (vol < 0)) { 
    vol = w;
  }

  findCM<t_scalar,brr3>(ips, W, cm);

  // Finally, get dr: 
  arr3Initialise<brr3>(cmA);
  arr3Initialise<brr3>(cmB);
  int cntA = 0, cntB = 0;
  for (int i=0; i<4; i++){
    arr3Resize<t_scalar,brr3>(countAin[i], tetAin[i]);
    arr3arr3Add<t_scalar,brr3>(tetAin[i], cmA, cmA);
    cntA += countAin[i];

    arr3Resize<t_scalar,brr3>(countBin[i], tetBin[i]);
    arr3arr3Add<t_scalar,brr3>(tetBin[i], cmB, cmB);
    cntB += countBin[i];
  } 
  if (cntA > 1) arr3Resize<t_scalar,brr3>(-1./cntA, cmA);
  if (cntB > 1) arr3Resize<t_scalar,brr3>(1./cntB, cmB);
  // - calculate dr:
  if (cntA > 0) {
    if (cntB == 0) arr3Store<t_scalar,brr3>(cmA, dr);
    else {
      arr3arr3Add<t_scalar,brr3>(cmA, cmB, dr); 
      arr3Resize<t_scalar,brr3>(0.5, dr); 
    }
  } else {
    arr3Store<t_scalar,brr3>(cmB,dr);
  }

  return vol; 
 

} 

/** Return the volume and area of intersection between two tetrahedra */ 
//  input: arr3 (&tetA)[4], arr3 (&tetB)[4])
//  output: scalar vol and scalar area. 
template <class t_scalar, class brr3> void volumeAndAreaIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4], t_scalar &vol, t_scalar &area){


  vol = 0.0; 
  t_scalar v_aux, a_aux;
  area = 0.0;
  brr3 W[56]; 

  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists<t_scalar,brr3>(tetA[i], ips, W) && nodeInTet<t_scalar,brr3>(tetA[i], tetB)) { 
      // cout << "tetA-node[" << i << "] found int tetB " << endl; 
      arr3Store<t_scalar,brr3>(tetA[i], W[ips]);
      ips += 1; 
      volumeAndAreaForNode<t_scalar,brr3>(tetA, i, v_aux, a_aux); 
      vol += v_aux;
      area += a_aux; 
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists<t_scalar,brr3>(tetB[i], ips, W) && nodeInTet<t_scalar,brr3>(tetB[i], tetA)) { 
      // cout << "tetB-node[" << i << "] found int tetA " << endl; 
      arr3Store<t_scalar,brr3>(tetB[i], W[ips]);
      ips += 1;
      volumeAndAreaForNode<t_scalar,brr3>(tetB, i, v_aux, a_aux); 
      vol += v_aux;
      area += a_aux; 
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  brr3 ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 0, 1, 2, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 1, 2, 3, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 2, 3, 0, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetA, i, j, tetB, 3, 0, 1, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 0, 1, 2, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 1, 2, 3, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 2, 3, 0, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint<t_scalar,brr3>(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists<t_scalar,brr3>(ip, ips, W)) {
          arr3Store<t_scalar,brr3>(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint<t_scalar,brr3>(ip, tetB, i, j, tetA, 3, 0, 1, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  area *= ffea_const::half; 
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) {
    vol = 0.;
    area = 0.;
    return; 
  } 
 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  maxVolumeAndArea<t_scalar,brr3>(ips, W, v_aux, a_aux);
  if ((fabs(vol) > v_aux) || (vol < 0)) { 
    vol = v_aux;
    area = a_aux;
  }

} 


// INSTANTIATE everything to arr3 and eventually to grr3: 
template bool exists<scalar,arr3>(arr3 &p, int ips, arr3 (&W)[56]);

template scalar maxVolume<scalar,arr3>(int ips, arr3 (&W)[56]);

template void maxVolumeAndArea<scalar,arr3>(int ips, arr3 (&W)[56], scalar &volume, scalar &area);

template void findCM<scalar,arr3>(int ips, arr3 (&W)[56], arr3 &cm);

template void getBAndN_Order<scalar,arr3>(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &t, arr3 &b, arr3 &n);

template void getBAndN<scalar,arr3>(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 &t, arr3 &b, arr3 &n);

template scalar volumeForNode<scalar,arr3>(arr3 (&tetA)[4], int node);

template void volumeAndAreaForNode<scalar,arr3>(arr3 (&tetA)[4], int node, scalar &volume, scalar &area);

template scalar volumeForIntPoint<scalar,arr3>(arr3 &ip, arr3 (&tetA)[4], int e1, int e2, arr3 (&tetB)[4], int f1, int f2, int f3);

template void volumeAndAreaForIntPoint<scalar,arr3>(arr3 &ip, arr3 (&tetA)[4], int e1, int e2, arr3 (&tetB)[4], int f1, int f2, int f3, scalar &volume, scalar &area);

template scalar volumeIntersection<scalar,arr3>(arr3 (&tetA)[4], arr3 (&tetB)[4]);
template scalar volumeIntersection<scalar,arr3>(arr3 (&tetA)[4], arr3 (&tetB)[4], arr3 &cm);
template scalar volumeIntersection<scalar,arr3>(arr3 (&tetA)[4], arr3 (&tetB)[4], arr3 &cm, arr3 &dr);

template void volumeAndAreaIntersection<scalar,arr3>(arr3 (&tetA)[4], arr3 (&tetB)[4], scalar &vol, scalar &area);

#ifndef USE_DOUBLE
template bool exists<geoscalar,grr3>(grr3 &p, int ips, grr3 (&W)[56]);
template geoscalar maxVolume<geoscalar,grr3>(int ips, grr3 (&W)[56]);
template void maxVolumeAndArea<geoscalar,grr3>(int ips, grr3 (&W)[56], geoscalar &volume, geoscalar &area);
template void findCM<geoscalar,grr3>(int ips, grr3 (&W)[56], grr3 &cm);
template void getBAndN_Order<geoscalar,grr3>(grr3 (&tetA)[4], int n0, int n1, int n2, grr3 &t, grr3 &b, grr3 &n);
template void getBAndN<geoscalar,grr3>(grr3 (&tetA)[4], int n0, int n1, int n2, grr3 &t, grr3 &b, grr3 &n);
template geoscalar volumeForNode<geoscalar,grr3>(grr3 (&tetA)[4], int node);
template void volumeAndAreaForNode<geoscalar,grr3>(grr3 (&tetA)[4], int node, geoscalar &volume, geoscalar &area);
template geoscalar volumeForIntPoint<geoscalar,grr3>(grr3 &ip, grr3 (&tetA)[4], int e1, int e2, grr3 (&tetB)[4], int f1, int f2, int f3);
template void volumeAndAreaForIntPoint<geoscalar,grr3>(grr3 &ip, grr3 (&tetA)[4], int e1, int e2, grr3 (&tetB)[4], int f1, int f2, int f3, geoscalar &volume, geoscalar &area);
template geoscalar volumeIntersection<geoscalar,grr3>(grr3 (&tetA)[4], grr3 (&tetB)[4]);
template geoscalar volumeIntersection<geoscalar,grr3>(grr3 (&tetA)[4], grr3 (&tetB)[4], grr3 &cm);
template geoscalar volumeIntersection<geoscalar,grr3>(grr3 (&tetA)[4], grr3 (&tetB)[4], grr3 &cm, grr3 &dr);
template void volumeAndAreaIntersection<geoscalar,grr3>(grr3 (&tetA)[4], grr3 (&tetB)[4], geoscalar &vol, geoscalar &area);
#endif 
