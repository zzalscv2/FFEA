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

#ifndef VECTORS_H
#define VECTORS_H

#include <cmath>

using namespace std;

class vector3
{
	public:

		// Constructors/Destructors
		vector3() {
			x = 0;
			y = 0;
			z = 0;
			r = 0;
		}

		vector3(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		}

		~vector3() {
			x = 0;
			y = 0;
			z = 0;
			r = 0;
		}

		// Overloaded operators
		vector3 operator+(vector3 rvalue) {
			vector3 result;
			result.x = x + rvalue.x;
			result.y = y + rvalue.y;
			result.z = z + rvalue.z;
			result.r = sqrt(pow(result.x, 2) + pow(result.y, 2) + pow(result.z, 2));
			return result;
		}

		void operator+=(vector3 rvalue) {
			x += rvalue.x;
			y += rvalue.y;
			z += rvalue.z;
			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		}

		vector3 operator-(vector3 rvalue) {
			vector3 result;
			result.x = x - rvalue.x;
			result.y = y - rvalue.y;
			result.z = z - rvalue.z;
			result.r = sqrt(pow(result.x, 2) + pow(result.y, 2) + pow(result.z, 2));
			return result;
		}

		void operator-=(vector3 rvalue) {
			x -= rvalue.x;
			y -= rvalue.y;
			z -= rvalue.z;
			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		}

		vector3 operator*(double rvalue) {
			vector3 result;
			result.x = x * rvalue;
			result.y = y * rvalue;
			result.z = z * rvalue;
			result.r = r * rvalue;
			return result;
		}

		void operator*=(double rvalue) {
			x *= rvalue;
			y *= rvalue;
			z *= rvalue;
			r *= rvalue;
		}

		// Member functions
		double get_mag() {
			return r;
		}

		double check_mag() {
			return r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		}
		
		void normalise() {
			x /= r;
			y /= r;
			z /= r;
			r = 1.0;
		}

		void set_pos(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		}
		
		double dot(vector3 rvec) {
			return x * rvec.x + y * rvec.y + z * rvec.z;
		}
		
		vector3 cross(vector3 rvec) {
			
			vector3 c;
			c.x = y * rvec.z - z * rvec.y;
			c.y = z * rvec.x - x * rvec.z;
			c.z = x * rvec.y - y * rvec.x;
			c.r = sqrt(pow(c.x, 2) + pow(c.y, 2) + pow(c.z, 2));
			return c;	
		}
		
		// Data members
		double x, y, z, r;
};
#endif
