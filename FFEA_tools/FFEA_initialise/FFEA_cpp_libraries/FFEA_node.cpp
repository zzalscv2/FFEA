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
