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

#ifndef FFEA_NODE_INCLUDED
#define FFEA_NODE_INCLUDED

#include <iostream>
#include "FFEA_math.hpp"

using namespace std;

class FFEA_node {

	public:

		// Constructor / Destructor
		FFEA_node();
		~FFEA_node();

		// Loading Functions
		int load(string fname);
		int load_from_node(string fname);
		int load_from_vol(string fname);

		// Data access functions
		int get_num_nodes();
		int get_num_interior_nodes();
		int get_num_surface_nodes();

		float get_centroid();

		// Member variables
		vector3f *pos;
		
	private:

		int num_nodes, num_interior_nodes, num_surface_nodes;
};
#endif
