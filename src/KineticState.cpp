#include "KineticState.h"

KineticState::KineticState() {
	conformation_index = 0;
	bound_sites.clear();
}

KineticState::~KineticState() {
	conformation_index = 0;
	bound_sites.clear();
}

int KineticState::init() {

	// Set conformation index
	conformation_index = 0;

	// No sites to bind
	bound_sites.clear();

	return FFEA_OK;
}

int KineticState::init(int conf_index, int *bound_site_types, int num_bsite_types, BindingSite *binding_sites, int total_num_binding_sites) {

	// Set conformation index
	conformation_index = conf_index;

	// Set some bound binding sites for this state
	for(int i = 0; i < num_bsite_types; ++i) {
		if(bound_site_types[i] == 1) {
			for(int j = 0; j < total_num_binding_sites; ++j) {
				bound_sites.insert(&binding_sites[j]);
			}
		}
	}

	return FFEA_OK;
}

/*void KineticState::print_details() {

	cout << "Kinetic Site:" << endl << endl;
	cout << "Conformation Index = " << conformation_index << endl;
	cout << "Active binding sites = ";
	for(set<int>::iterator it = bound_site.begin(); it != bound_site.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
}*/

