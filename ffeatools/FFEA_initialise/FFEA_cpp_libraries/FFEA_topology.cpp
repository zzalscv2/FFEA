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

#include "FFEA_topology.hpp"


// Topology stuff
FFEA_topology::FFEA_topology() {
	
	num_elements = 0;
	element = NULL;
}

FFEA_topology::~FFEA_topology() {
	
	num_elements = 0;

	delete[] element;
	element = NULL;
}

int FFEA_topology::load(string fname) {

	// This function checks the extension, then sends it to the filetype specific functions
	
	return 0;
}

int FFEA_topology::load_from_top(string fname) {

	return 0;
}

int FFEA_topology::load_from_vol(string fname) {

	return 0;
}

int FFEA_topology::get_num_elements() {
	
	return num_elements;
}

// Individual element stuff
FFEA_element::FFEA_element() {

	eltype = LINEAR;
	index = NULL;
	pos = NULL;
}

FFEA_element::~FFEA_element() {

	eltype = LINEAR;
	
	delete[] index;
	index = NULL;

	delete[] pos;
	pos = NULL;
}

int FFEA_linear_element::set_indices(int *n) {

	if (index == NULL) {
		index = new int[4];
	}

	for(int i = 0; i < 4; ++i) {
		index[i] = n[i];
	}

	eltype = LINEAR;
}
