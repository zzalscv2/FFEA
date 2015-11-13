#ifndef KINETICBINDINGSITE_H_INCLUDED
#define KINETICBINDINGSITE_H_INCLUDED

#include <vector>
#include <set>
#include <iostream>
#include "Face.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"

using namespace std;

class KineticBindingSite {

	public:

		KineticBindingSite();
		~KineticBindingSite();

		int init(int site_type, int blob_index, int conf_index, vector<Face*> face_vector);
		void calc_num_nodes();
		void calc_centroid();
		void calc_area();

		int site_type;
		int blob_index;
		int conf_index;
		vector<Face*> faces;

		int num_nodes;
		vector3 centroid;
		scalar area;
};

#endif
