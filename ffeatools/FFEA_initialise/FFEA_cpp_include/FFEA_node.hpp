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
