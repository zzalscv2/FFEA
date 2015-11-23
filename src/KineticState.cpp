#include "KineticState.h"

KineticState::KineticState() {
	conformation_index = 0;
	active_site.clear();
}

KineticState::~KineticState() {
	conformation_index = 0;
	active_site.clear();
}

int KineticState::init(int conf_index, int *active_bsites, int num_bsite_types) {

	// Set conformation index
	conformation_index = conf_index;

	// Set some active binding sites
	for(int i = 0; i < num_bsite_types; ++i) {
		if(active_bsites[i] == 1) {
			active_site.push_back(i);
		}
	}

	return FFEA_OK;
}

void KineticState::print_details() {

	cout << "Kinetic Site:" << endl << endl;
	cout << "Conformation Index = " << conformation_index << endl;
	cout << "Active binding sites = ";
	for(vector<int>::iterator it = active_site.begin(); it != active_site.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
}

