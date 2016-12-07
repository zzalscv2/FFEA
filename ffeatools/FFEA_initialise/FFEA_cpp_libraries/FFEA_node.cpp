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

#include "FFEA_node.hpp"

// Constructor / Destructor
FFEA_node::FFEA_node() {
	num_nodes = 0;
	num_interior_nodes = 0;
	num_surface_nodes = 0;
	pos = NULL;
}

FFEA_node::~FFEA_node() {
	num_nodes = 0;
	num_interior_nodes = 0;
	num_surface_nodes = 0;
	pos = NULL;
}

// Loading functions
int FFEA_node::load(string fname) {

	// This function checks the extension, then sends it to the filetype specific functions
	
	return 0;
}

int FFEA_node::load_from_node(string fname) {

	return 0;
}

int FFEA_node::load_from_vol(string fname) {

	return 0;
}

// Data access functions
int FFEA_node::get_num_nodes() {

	return num_nodes;
}

int FFEA_node::get_num_interior_nodes() {

	return num_interior_nodes;
}

int FFEA_node::get_num_surface_nodes() {

	return num_surface_nodes;
}

float FFEA_node::get_centroid() {

	
}
