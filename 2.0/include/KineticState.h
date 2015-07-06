#ifndef KINETICSTATES_H_INCLUDED
#define KINETICSTATES_H_INCLUDED

#define FFEA_KINETIC_STATE_UNBOUND	0
#define FFEA_KINETIC_STATE_BOUND	1

#include "FFEA_return_codes.h"
#include <stdio.h>
#include <iostream>

using namespace std;

class KineticState {
	
	public:
		
		KineticState();

		~KineticState();

		int init(int conf_ind, int bound_state, int site_type_from, int site_type_to);

		// Active Conformation
		int conformation_index;

		// Bound or Unbound
		int bound;
		int binding_site_type_from;
		int binding_site_type_to;
		
};
#endif
