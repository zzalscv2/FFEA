#ifndef BINDINGMATRIXHINCLUDED
#define BINDINGMATRIXHINCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "Face.h"

using namespace std;

class BindingSite{

	public:

		BindingSite();
		~BindingSite();

		int num_faces, site_type;
		vector3 centroid;
		scalar area, length;
		vector<Face*> faces;
		
		void set_num_faces(int num_faces);
		void set_type(int site_type);
		void add_face(Face *aface);
		vector3 * calc_centroid();
		void calc_site_shape();
		scalar calc_size();
		scalar calc_area();
};

class BindingSite_matrix{

	public:
		
		BindingSite_matrix();
		~BindingSite_matrix();

		int init(string fname);
		int get_num_types();

		// Public variables
		int num_interaction_types;
		int **interaction;
};

#endif
