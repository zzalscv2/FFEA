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

#ifndef FFEA_TOPOLOGY_INCLUDED
#define FFEA_TOPOLOGY_INCLUDED


#include <iostream>
using namespace std;

// Define an enum for element types
enum elementType {LINEAR, SECONDARY};

class FFEA_element {

	public:

		// Constructor / Destructor
		FFEA_element();
		~FFEA_element();
		
	protected:

		int *index;
		float **pos;
		elementType eltype;
};


// We don't need constructors! Only the assignment functions need overloading...
class FFEA_linear_element: public FFEA_element {

	public:

		// Assignment functions
		int set_indices(int *n);
		int set_indices(int n0, int n1, int n2, int n3);

		int set_pos(float **pos);
		//int set_pos(float *pos0, *pos1, *pos2, *pos3);
};


class FFEA_secondary_element: public FFEA_element {

	public:

		// Assignment functions
		int set_indices(int *n);
		int set_indices(int n0, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9);

		int set_pos(float **pos);
		//int set_pos(float *pos0, *pos1, *pos2, *pos3, *pos, *pos5, *pos6, *pos7, *pos8, *pos9);	
};

class FFEA_topology {

	public:

		// Constructor / Destructor
		FFEA_topology();
		~FFEA_topology();

		// Loading Functions
		int load(string fname);
		int load_from_top(string fname);
		int load_from_vol(string fname);

		// Data access functions
		int get_num_elements();

		FFEA_element *element;

	private:

		int num_elements;
};

#endif
