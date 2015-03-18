#include "mat_vec_types.h"

class Spring {

	public:
	
		Spring() {
			blob_index = new int[2];
			conformation_index = new int[2];
			node_index = new int[2];
			k = 0.0;
			l = 0.0;
			am_i_active = false;
		}

		~Spring() {
			delete[] blob_index;
			blob_index = NULL;
			delete[] conformation_index;
			conformation_index = NULL;
			delete[] node_index;
			node_index = NULL;
			k = 0.0;
			l = 0.0;
			am_i_active = false;
		}

		/* 
		 *  Variables
		 */

		// Spring constant
		scalar k;
		
		// Equilibrium length
		scalar l;
		
		// Blobs connected to
		int *blob_index;

		// Conformations connected to
		int *conformation_index;
	
		// Nodes connected to
		int *node_index;

		// Check if spring is active
		bool am_i_active;
};
