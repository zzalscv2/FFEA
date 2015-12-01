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

		// Structure
		int num_faces;
		vector<Face*> faces;

		// Properties
		vector3 centroid;
		scalar area, length;
		
		void set_num_faces(int num_faces);
		void set_type(int site_type);
		int get_type();
		void add_face(Face *aface);
		void calc_centroid(vector3 cent);
		void calc_site_shape(vector3 cent);
		scalar calc_size();
		scalar calc_area();

	private:

		int site_type;
};

class BindingSite_matrix{

	public:
		
		BindingSite_matrix();
		~BindingSite_matrix();

		int init(string fname);
		int get_num_interaction_types();
		void print_to_screen();

		// Public variables
		int num_interaction_types;
		int **interaction;
};

#endif
