#include "KineticState.h"

KineticState::KineticState() {
	conformation_index = 0;
	bound = FFEA_KINETIC_STATE_UNBOUND;
	binding_site_type_from = -1;
	binding_site_type_to = -1;
}

KineticState::~KineticState() {
	conformation_index = 0;
	bound = FFEA_KINETIC_STATE_UNBOUND;
	binding_site_type_from = -1;
	binding_site_type_to = -1;
}

int KineticState::init(int conf_ind, int bound_state, int site_type_from, int site_type_to) {

	conformation_index = conf_ind;

	if(bound_state != FFEA_KINETIC_STATE_UNBOUND && bound_state != FFEA_KINETIC_STATE_BOUND) {
		FFEA_error_text();
		cout << "In 'KineticState::init', expected bound_state = 0 or bound_state = 1. Got bound_state = " << bound_state << "." << endl;
		return FFEA_ERROR;
	} else {
		bound = bound_state;
	}

	binding_site_type_from = site_type_from;
	binding_site_type_to = site_type_to;
	return FFEA_OK;
}

