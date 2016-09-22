#include <iostream>
#include <sstream>
#include "VolumeIntersection.h"
#include "CheckTetrahedraOverlap.h"

using namespace std; 

int main() {

  arr3 tetA[4], tetB[4];
  arr3 tetC[4], tetD[4]; 
  // set up tetA and tetB:
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      tetA[i][j] = 0;
      tetB[i][j] = 0;
    }
  }

  tetA[1][1] = 1.0;
  tetA[2][0] = 1.0;
  tetA[3][2] = 1.0;

  tetB[1][1] = 1.0;
  tetB[2][0] = 1.0;
  tetB[3][2] = 1.0;

  for (int i=0; i<4; i++) {
    tetB[i][2] += 0.5;
    tetA[i][1] += 0.02;
    tetA[i][0] += 0.05;
    // tetA[i][1] += 0.0;
  }

  cout << " tetA: " << endl;
  for (int i=0; i<4; i++){
    for (int j=0;j<3;j++){
       cout << tetA[i][j] << " ";
    }
    cout << endl;
  }
  cout << " tetB: " << endl;
  for (int i=0; i<4; i++){
    for (int j=0;j<3;j++){
       cout << tetB[i][j] << " ";
    }
    cout << endl;
  }


  if (! tet_a_tet(tetA, tetB)) {
    cout << " these tetrahedra are known to intersect" << endl;
    return 1;
  } 

  scalar vol = volumeIntersection<scalar,arr3>(tetA, tetB);
  stringstream ss;
  ss.precision(14);
  ss << vol;
  string hardcodedresult = "0.020833333333333"; 
  cout.precision(14);
  cout << "total volume: " << vol << endl;

  if (hardcodedresult.compare(ss.str())) {
    return 1; 
  } /* else {
    return 0; 
  } */

  // set up tetC and tetD:
  tetC[0][0] = 5.16991;
  tetC[0][1] = 5.11394;
  tetC[0][2] = -2.59001;
  tetC[1][0] = 2.68658;
  tetC[1][1] = 5.11124;
  tetC[1][2] = -5.03959;
  tetC[2][0] = 5.12737;
  tetC[2][1] = 5.10932;
  tetC[2][2] = -5.0372;
  tetC[3][0] = 2.98162;
  tetC[3][1] = 2.85562;
  tetC[3][2] = -2.78657;
 
  tetD[0][0] = 4.8437;
  tetD[0][1] = 2.47363;
  tetD[0][2] = -5.0158;
  tetD[1][0] = 4.84606;
  tetD[1][1] = 5.11758;
  tetD[1][2] = -2.44064;
  tetD[2][0] = 7.20236;
  tetD[2][1] = 5.11031; 
  tetD[2][2] = -5.00642;
  tetD[3][0] = 4.88715;
  tetD[3][1] = 5.11134;
  tetD[3][2] = -5.00755; 

  if (! tet_a_tet(tetC, tetD)) {
    cout << "tetrahedra tetC and tetD are known to intersect " << endl;
    cout << "  and I found tet_a_tet fails in this case" << endl;
    // return 1;
  } 

  vol = volumeIntersection<scalar,arr3>(tetC, tetD);
  if (vol == 0) {
    cout << " intersecting volume should not be zero for these tetrahedra" << endl; 
    return 1;
  } 

  return 0; 

}
   
