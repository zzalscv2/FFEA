#ifndef KINETICSTATES_H_INCLUDED
#define KINETICSTATES_H_INCLUDED

#define FFEA_KINETIC_STATE_UNBOUND	0
#define FFEA_KINETIC_STATE_BOUND	1

#include "FFEA_return_codes.h"
#include "BindingSite.h"

#include <stdio.h>
#include <iostream>
#include <set>
using namespace std;

class KineticState {
	
	public:
		
		KineticState();

		~KineticState();

		int init();
		int init(int conf_index, int *bound_site_types, int num_bsite_types, BindingSite *binding_sites, int total_num_binding_sites);

		//void print_details();

		// Active Conformation
		int conformation_index;

		// Set of bound binding sites
		set<BindingSite*> bound_sites;
};
#endif
