#include "BindingSite.h"

BindingSite_matrix::BindingSite_matrix() {

	num_interaction_types = 0;
	interaction = NULL;
}

BindingSite_matrix::~BindingSite_matrix() {

	num_interaction_types = 0;
	interaction = NULL;
}

int BindingSite_matrix::init(string fname) {

	// Open file
	ifstream fin;
	fin.open(fname.c_str());
	if(fin.fail()) {
		FFEA_ERROR_MESSG("'binding_params_fname' %s not found\n", fname.c_str())
	}

	// Check if correct file
	int MAX_BUF_SIZE = 255;
	char buf[MAX_BUF_SIZE];
	string buf_string;
	fin.getline(buf, MAX_BUF_SIZE);
	buf_string = string(buf);
	boost::trim(buf_string);
	if(buf_string != "ffea binding site params file") {
		FFEA_ERROR_MESSG("Expected 'ffea binding site params file', got '%s'\n", buf_string.c_str())
	}

	// Get num_site_types
	fin >> buf_string >> num_interaction_types;
	
	// Get all interactions
	interaction = new int*[num_interaction_types];
	for(int i = 0; i < num_interaction_types; ++i) {
		interaction[i] = new int[num_interaction_types];
		for(int j = 0; j < num_interaction_types; ++j) {
			fin >> interaction[i][j];
		}
	}
	fin.close();
	return FFEA_OK;
}

int BindingSite_matrix::get_num_types() {
	
	return num_interaction_types;
}

BindingSite::BindingSite() {
	num_faces = 0;
	site_type = -1;
	faces.clear();
	vector3_set_zero(&centroid);
	area = 0.0;
	length = 0.0;
}

BindingSite::~BindingSite() {
	num_faces = 0;
	site_type = -1;
	faces.clear();
	vector3_set_zero(&centroid);
	area = 0.0;
	length = 0.0;
}

void BindingSite::set_num_faces(int num_faces) {
	this->num_faces = num_faces;
}

void BindingSite::set_type(int site_type) {
	this->site_type = site_type;
}

void BindingSite::add_face(Face *aface) {
	faces.push_back(aface);
}

vector3 * BindingSite::calc_centroid() {
	
	vector3_set_zero(&centroid);
	vector3 *face_centroid;
	vector<Face*>::iterator it;
	for(it = faces.begin(); it != faces.end(); ++it) {
		face_centroid = (*it)->get_centroid();
		centroid.x += face_centroid->x;
		centroid.y += face_centroid->y;
		centroid.z += face_centroid->z;
	}
	vec3_scale(&centroid, 1.0/num_faces);
	return &centroid;
}

scalar BindingSite::calc_area() {

	return area;
}

scalar BindingSite::calc_size() {
	
	calc_area();
	length = sqrt(area / M_PI);
	return length;
}

void BindingSite::calc_site_shape() {
	
	calc_centroid();
	calc_size();
}