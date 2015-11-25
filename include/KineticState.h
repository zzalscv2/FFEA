#ifndef KINETICSTATES_H_INCLUDED
#define KINETICSTATES_H_INCLUDED

#define FFEA_KINETIC_STATE_UNBOUND	0
#define FFEA_KINETIC_STATE_BOUND	1

#include "FFEA_return_codes.h"
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

class KineticState {
	
	public:
		
		KineticState();

		~KineticState();

		int init(int conf_index, int *active_bsites, int num_bsite_types);

		void print_details();

		// Active Conformation
		int conformation_index;

		// Set of active binding site types
		vector<int> active_site;

		// And how many are active?
		int num_active_bsites;	
};
#endif
