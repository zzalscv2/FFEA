#include "KineticState.h"

KineticState::KineticState() {
	conformation_index = 0;
	base = -1;
	target = -1;
	bound = false;
	base_site = NULL;
	target_site = NULL;
}

KineticState::~KineticState() {
	conformation_index = 0;
	base = -1;
	target = -1;
	bound = false;
	base_site = NULL;
	target_site = NULL;
}

int KineticState::init() {

	// Set conformation index
	conformation_index = 0;

	// No sites to bind
	base = -1;
	target = -1;
	bound = false;

	return FFEA_OK;
}

int KineticState::init(int conf_index, int base, int target) {

	// Set conformation index
	conformation_index = conf_index;

	// Set the binding site types (from type 'base' to type 'target')
	this->base = base;
	this->target = target;

	if(base == -1 && target == -1) {
		bound = false;
	} else if (base == -1 && target != -1 || base != -1 && target == -1) {
		FFEA_ERROR_MESSG("Either both site types must be -1 (not bound) or neither should be -1 (bound)")
	} else {
		bound = true;
	}

	return FFEA_OK;
}

int KineticState::get_conformation_index() {
	
	return conformation_index;
}

bool KineticState::is_bound() {

	return bound;
}

int KineticState::get_base_bsite_type() {

	return base;
}

int KineticState::get_target_bsite_type() {

	return target;
}

void KineticState::set_sites(BindingSite *base, BindingSite *target) {

	base_site = base;
	target_site = target;
}

BindingSite * KineticState::get_base_site() {
	
	return base_site;
}

BindingSite * KineticState::get_target_site() {
	
	return target_site;
}
