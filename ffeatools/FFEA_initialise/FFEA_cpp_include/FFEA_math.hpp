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

#ifndef FFEA_MATH_INCLUDED
#define FFEA_MATH_INCLUDED

#include <iostream>
using namespace std;

class vector3i {

	public:

		// Constructor / Destructor
		vector3i();
		~vector3i();

		// Overload Operators
		vector3i operator+(vector3i rvalue);
		vector3i operator+=(vector3i rvalue);
		vector3i operator-(vector3i rvalue);
		vector3i operator-=(vector3i rvalue);
		vector3i operator*(float rvalue);
		vector3i operator*=(float rvalue);

		// Data access functions
		void set_components(int a, int b, int c);
		float get_magnitude();

		// Operation functions
		void normalise();

		// Member variables
		int comp[3];
		
	private:

		float magnitude;
};

class vector3f {

	public:

		// Constructor / Destructor
		vector3f();
		~vector3f();

		// Overload Operators
		vector3f operator+(vector3f rvalue);
		vector3f operator+=(vector3f rvalue);
		vector3f operator-(vector3f rvalue);
		vector3f operator-=(vector3f rvalue);
		vector3f operator*(float rvalue);
		vector3f operator*=(float rvalue);

		// Data access functions
		void set_components(float a, float b, float c);
		float get_magnitude();

		// Operation functions
		void normalise();

		// Member variables
		float comp[3];
		
	private:

		float magnitude;
};
#endif
