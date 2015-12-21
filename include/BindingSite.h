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
		int num_faces, site_type;

		/** List of faces making up the site */
		vector<Face*> faces;
		
		void print_to_screen();
		void set_num_faces(int num_faces);
		void set_type(int site_type);
		int get_type();
		void add_face(Face *aface);
		set<int> get_nodes();

		void calculate_centroid();
		vector3 get_centroid();
		void calculate_area();
		scalar get_area();
		void calculate_characteristic_length();
		scalar get_characteristic_length();

		static bool sites_in_range(BindingSite a, BindingSite b);

	private:

		/** Centroid of the whole site (needs recalculating if simulation has continued) */
		vector3 centroid;

		/** Area of the whole site (needs recalculating if simulation has continued) */
		scalar area;

		/** Area of the whole site (needs recalculating if area has updated */
		scalar characteristic_length;
};

class BindingSite_matrix{

	public:
		
		BindingSite_matrix();
		~BindingSite_matrix();

		int init(string fname);
		int get_num_interaction_types();
		void print_to_screen();

		/** Number of interaction types allowed */
		int num_interaction_types;
		
		/** 2D matrix defining allowed interactions between site types */
		bool **interaction;
};

#endif
