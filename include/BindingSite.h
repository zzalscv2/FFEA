#ifndef BINDINGMATRIXHINCLUDED
#define BINDINGMATRIXHINCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
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
		scalar area, radius;
		
		void set_num_faces(int num_faces);
		void set_type(int site_type);
		int get_type();
		void add_face(Face *aface);
		vector3 calc_centroid();
		void calc_dimensions();
		scalar calc_size();
		scalar calc_area();
		set<int> get_nodes();
		vector3 get_centroid();

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
