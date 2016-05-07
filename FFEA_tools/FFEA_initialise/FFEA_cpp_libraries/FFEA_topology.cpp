#include "FFEA_topology.h"


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

	eltype = LINEAR_TYPE
	index = NULL;
	pos = NULL;
}

FFEA_element::~FFEA_element() {

	eltype = LINEAR_TYPE
	
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
		try {
			index = 
		} catch {

		}
	}
}
