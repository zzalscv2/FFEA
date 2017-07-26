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

typedef double scalar;
typedef scalar arr3[3];


class vector3 {
public:
    arr3 data;
    scalar& x = data[0];
    scalar& y = data[1];
    scalar& z = data[2];
    scalar& operator [](std::size_t i) { return data[i]; }
    void assign( scalar t0, scalar t1, scalar t2 )
    { data[0] = t0; data[1] = t1; data[2] = t2; }
};

int main() {

  vector3 v;
  v.assign(0, 1, 2);
  std::cout << v[0] << ", " << v[1] << ", " << v[2] << std::endl; 

  return 0;

}
