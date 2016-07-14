#include <iostream>
#include <sstream>
#include "VolumeIntersection.h"
#include "CheckTetrahedraOverlap.h"

using namespace std; 

int main() {

  arr3 tetA[4], tetB[4];
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
  } else {
    return 0; 
  } 

}
   
