#include "KineticBindingSite.h"

KineticBindingSite::KineticBindingSite() {
	site_type = 0;
	blob_index = 0;
	conf_index = 0;
	faces.clear();

	num_nodes = 0;
	centroid.x = 0.0;
	centroid.y = 0.0;
	centroid.z = 0.0;
	area = 0.0;
}

KineticBindingSite::~KineticBindingSite() {
	site_type = 0;
	blob_index = 0;
	conf_index = 0;
	faces.clear();

	num_nodes = 0;
	centroid.x = 0.0;
	centroid.y = 0.0;
	centroid.z = 0.0;
	area = 0.0;
}

int KineticBindingSite::init(int site_type, int blob_index, int conf_index, vector<Face*> face_vector) {
	
	// Init stuff
	vector<Face*>::iterator it;	
	this->site_type = site_type;
	this->blob_index = blob_index;
	this->conf_index = conf_index;
	for(it = face_vector.begin(); it != face_vector.end(); ++it) {
		faces.push_back(*it);
	}
	
	// Calculate stuff
	calc_num_nodes();
	calc_centroid();
	calc_area();
	return FFEA_OK;
}

void KineticBindingSite::calc_num_nodes() {
	
	int i;
	vector<Face*>::iterator it;
	set<int> node_set;
	for(it = faces.begin(); it != faces.end(); ++it) {
		for(i = 0; i < 3; ++i) {
			node_set.insert((*it)->n[i]->index);
		}
	}
	num_nodes = node_set.size();
}

void KineticBindingSite::calc_centroid() {

	vector<Face*>::iterator it;
	centroid.x = 0.0;
	centroid.y = 0.0;
	centroid.z = 0.0;

	vector3 *face_centroid;
	for(it = faces.begin(); it != faces.end(); ++it) {
		face_centroid = (*it)->get_centroid();
		centroid.x += 3 * face_centroid->x;
		centroid.y += 3 * face_centroid->y;
		centroid.z += 3 * face_centroid->z;
	}

	centroid.x *= 1.0 / num_nodes;
	centroid.y *= 1.0 / num_nodes;
	centroid.z *= 1.0 / num_nodes;

}

void KineticBindingSite::calc_area() {

	vector<Face*>::iterator it;
	area = 0.0;
	for(it = faces.begin(); it != faces.end(); ++it) {
		
	}
}

int KineticBindingSite::get_type() {

	return site_type;
}
