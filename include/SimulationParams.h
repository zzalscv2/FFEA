#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"

#define WALL_TYPE_PBC 0
#define WALL_TYPE_HARD 1
#define WALL_TYPE_STOP 2

#define MAX_FNAME_SIZE 200

/*
 * Simulation parameters
 */
class SimulationParams
{
	public:
		scalar dt;		/* time step */
		long long num_steps;		/* Number of time steps to run simulation for */
		int check;		/* Every how many steps should the program 'check' the system
					   i.e calculate energies, print snapshots etc. */
		int num_blobs;		/* Number of blobs in the system */
		int *num_conformations;	/* Number of conformations for each blob */
		int *num_states;	/* Number of states for each blob */
		int rng_seed;		/* Seed for random number generator */

		scalar kT;		/* boltzmann's const times temperature */

		int max_iterations_cg;	/* Max number of iterations when using conjugate gradient solver */
		scalar epsilon2;	/* The tolerance threshold for CG solver (squared) */

		int es_update;		/* Every how many steps should the electrostatic potential be recalculated */
		int es_N_x;		/* X dimension of the 3D lookup grid (in number of cells) */
		int es_N_y;		/* Y dimension of the 3D lookup grid (in number of cells) */
		int es_N_z;		/* Z dimension of the 3D lookup grid (in number of cells) */
		scalar es_h;		/* Dimension of each cell in the lookup grid (in multiples of inverse kappa) */

		scalar kappa;		/* Inverse Debye Screening length */

		scalar epsilon_0;	/* Permittivity of free space */
		scalar dielec_ext;	/* Exterior dielectric constant */

		int restart;		/* Whether or not to restart the simulation from the last available time step */

		int calc_vdw;		/* Whether or not to simulate van der waals interactions between surfaces */
		int calc_es;		/* Whether or not to simulate electrostatic interactions between proteins */
		int calc_noise;  	/* Whether or noise to simulate thermal noise for the system. I would, but it's up to you */

		int wall_x_1;
		int wall_x_2;
		int wall_y_1;
		int wall_y_2;
		int wall_z_1;
		int wall_z_2;

		int sticky_wall_xz;

		int do_stokes;
		scalar stokes_visc;

		scalar vdw_r_eq, vdw_eps;

		char trajectory_out_fname[MAX_FNAME_SIZE];
		char **measurement_out_fname;
		char temp_fname[MAX_FNAME_SIZE];
		char stress_out_fname[MAX_FNAME_SIZE];
		int stress_out_fname_set;

		char vdw_params_fname[MAX_FNAME_SIZE];

		SimulationParams();

		~SimulationParams();

		int validate();

		/*
		 * Expects a string of the form "lvalue = rvalue" where lvalue must be a recognised parameter name.
		 */
		int parse_param_assignment(char *str);

	private:
		int trajectory_out_fname_set;
		int measurement_out_fname_set;
		int vdw_params_fname_set;
};

#endif
