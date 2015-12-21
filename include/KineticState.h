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
		int init(int conf_index, int from, int to);

		int get_conformation_index();
		int get_base_bsite_type();
		int get_target_bsite_type();
		bool is_bound();
		void set_sites(BindingSite *base, BindingSite *target);
		BindingSite * get_base_site();
		BindingSite * get_target_site();

	private:

		/** @brief The active conformation in this state */
		int conformation_index;

		/** @brief The active binding site type in this blob in this state */
		int base;

		/** @brief The target binding site type on any blob in this state */
		int target;

		/** @brief Is the base site type bound in this state? */
		bool bound;

		/** @brief Pointers to the actual binding sites (after state change has occured) */
		BindingSite *base_site, *target_site;
};
#endif
