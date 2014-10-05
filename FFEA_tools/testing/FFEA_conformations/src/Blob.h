/*
 *      Blob.h
 */

#ifndef BLOB_H_INCLUDED
#define BLOB_H_INCLUDED

#include "SparsityPattern.h"
#include "ConnectivityTypes.h"

#include <stdio.h>
#include <math.h>

#include <omp.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "ConjugateGradientSolver.h"
#include "MassLumpedSolver.h"
#include "NoMassCGSolver.h"
#include "Face.h"
#include "CG_solver.h"
#include "Blob.h"
#include "BEM_Poisson_Boltzmann.h"
#include "LJ_matrix.h"

/*
 * The "Blob" class
 */
class Blob {

	public:
		/*
		 * Blob constructor:
		 * Initialises all variables and pointers to 0 (or NULL). Does not perform any memory allocation. Actual Blob initialisation
		 * is carried out by the init() method.
		 */
		Blob()
		{
			/* Initialise everything to zero */
			num_nodes = 0;
			num_elements = 0;
			num_surface_faces = 0;
			num_surface_nodes = 0;
			num_interior_nodes = 0;
			num_surface_elements = 0;
			num_interior_elements = 0;
			mass = 0;
			vdw_bb_energy = 0.0;
			blob_state = FFEA_BLOB_IS_STATIC;
			node = NULL;
			elem = NULL;
			surface = NULL;
			solver = NULL;
			linear_solver = 0;
			force = NULL;
			rng = NULL;
			poisson_solver = NULL;
			phi_Omega = NULL;
			phi_Gamma = NULL;
			q = NULL;
			nodal_q = NULL;
			poisson_rhs = NULL;
			num_pinned_nodes = 0;
			pinned_nodes_list = NULL;
		}

		/*
 		* Blob destructor:
 		* Closes all open files, deallocates memory held by arrays, solvers etc. and sets everything to zero (NULL)
 		*/
		~Blob() {
			/* Release the node, element and surface arrays */
			delete[] node;
			node = NULL;
			delete[] elem;
			elem = NULL;
			delete[] surface;
			surface = NULL;

			/* Release the force vector */
			delete[] force;
			force = NULL;

			/* Release the Solver */
			delete solver;
			solver = NULL;
			linear_solver = 0;

			/* Release the Poisson Solver */
			delete poisson_solver;
			poisson_solver = NULL;
			delete[] phi_Omega;
			phi_Omega = NULL;
			delete[] phi_Gamma;
			phi_Gamma = NULL;
			delete[] q;
			q = NULL;
			delete[] nodal_q;
			nodal_q = NULL;
			delete[] poisson_rhs;
			poisson_rhs = NULL;
			delete[] pinned_nodes_list;
			pinned_nodes_list = NULL;

			/* Set relevant data to zero */
			num_nodes = 0;
			num_elements = 0;
			num_surface_elements = 0;
			num_interior_elements = 0;
			num_surface_faces = 0;
			num_surface_nodes = 0;
			num_interior_nodes = 0;
			blob_state = FFEA_BLOB_IS_STATIC;
			mass = 0;
			vdw_bb_energy = 0.0;
			rng = NULL;
			num_pinned_nodes = 0;
		}

		/*
		 * Allocates the node, element and embedded charge arrays, initialised with the data
		 * read from the given node, element and embedded charge files.
		 * 'linear_solver' sets which type of linear solver is to be used: 0 for direct (forward/backward
		 * substitution) and 1 for iterative (preconditioned gonjugate gradient).
		 * Also takes the simulation parameters and the array of RNGs (for multiprocessor runs).
		 */
		int init( const char *node_filename, const char *topology_filename, const char *surface_filename, const char *material_params_filename,
				const char *stokes_filename, const char *vdw_filename, const char *pin_filename, scalar scale, int linear_solver,
				int blob_state, SimulationParams *params, LJ_matrix *lj_matrix, MTRand rng[], int num_threads)
		{

			//Load the node, topology, surface, materials and stokes parameter files.
			if(load_nodes(node_filename, scale) == FFEA_ERROR) {
				FFEA_ERROR_MESSG("Error when loading Blob nodes.\n");
			}
			if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
				if(load_topology(topology_filename) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob topology.\n");
				}
				if(load_surface(surface_filename) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob surface.\n");
				}
			} else {
				if(load_surface_no_topology(surface_filename) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob surface (ignoring topology).\n");
				}
			}
			if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
				if(load_material_params(material_params_filename) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob element material params.\n");
				}
				if(load_stokes_params(stokes_filename, scale) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob stokes parameter file.\n");
				}
			}

			this->lj_matrix = lj_matrix;

			if(load_vdw(vdw_filename, lj_matrix->get_num_types()) == FFEA_ERROR) {
				FFEA_ERROR_MESSG("Error when loading VdW parameter file.\n")
			}

			if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
				if(load_pinned_nodes(pin_filename) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when loading Blob pinned nodes file.\n");
				}
			} else {
				num_pinned_nodes = 0;
				pinned_nodes_list = NULL;
			}

			// Store the pointer to the simulation parameters class
			this->params = params;

			// Store the pointer to the random number generator array
			this->rng = rng;

			// Get the blob state
			this->blob_state = blob_state;

			// Need to know solver type
			this->linear_solver = linear_solver;

			// Linearise all the elements
			for(int i = 0; i < num_elements; i++) {
				elem[i].linearise_element();
			}

			// Get the rest jacobian, rest volume etc. of this Blob and store it for later use
			calc_rest_state_info();

			// Calculate the connectivity of the mesh, giving each node a list of pointers to elements
			// it is a member of. The pointers will be to the exact memory location in that element's struct
			// in which can be found the contribution to the force on that node.
			if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
				printf("\t\tCalculating node-element connectivity...");
				calculate_node_element_connectivity();
				printf("\t\tdone\n");

				// Create the chosen linear equation Solver for this Blob
				if(linear_solver == FFEA_DIRECT_SOLVER) {
					solver = new SparseSubstitutionSolver();
				} else if(linear_solver == FFEA_ITERATIVE_SOLVER) {
					solver = new ConjugateGradientSolver();
				} else if(linear_solver == FFEA_MASSLUMPED_SOLVER) {
					solver = new MassLumpedSolver();
				} else if(linear_solver == FFEA_NOMASS_CG_SOLVER) {
					solver = new NoMassCGSolver();
				} else {
					FFEA_ERROR_MESSG("Error in Blob initialisation: linear_solver=%d is not a valid solver choice\n", linear_solver);
				}

				// Initialise the Solver (whatever it may be)
				printf("\t\tBuilding solver:\n");
				if(solver->init(num_nodes, num_elements, node, elem, params, num_pinned_nodes, pinned_nodes_list) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error initialising solver.\n")
				}
			}


			// Allocate the force vector array for the whole Blob
			force = new vector3[num_nodes];
        		for(int i = 0; i < num_nodes; i++) {
                		force[i].x = 0;
                		force[i].y = 0;
                		force[i].z = 0;
        		}

			// Calculate how many faces each surface node is a part of
			num_contributing_faces = new int[num_surface_nodes];
			for(int i = 0; i < num_surface_nodes; i++) {
				num_contributing_faces[i] = 0;
			}
			for(int i = 0; i < num_surface_faces; i++) {
				for(int j = 0; j < 3; j++) {
					num_contributing_faces[surface[i].n[j]->index]++;
				}
			}

			// Store stokes drag on nodes, for use in viscosity matrix
			if(params->do_stokes == 1) {
				for(int i = 0; i < num_nodes; ++i) {
					node[i].stokes_drag = 6.0 * 3.141592654 * params->stokes_visc * node[i].stokes_radius;
				}
			}

			// Generate the Mass matrix for this Blob
			if(params->calc_es == 1) {
				build_mass_matrix();
			}

			// Create and initialise the poisson solver
			if(num_interior_nodes > 0) {

				// Calculate the Sparsity Pattern for the Poisson matrix (interior and exterior)
				printf("\t\tCalculating sparsity pattern for Poisson matrix and RHS 'knowns' matrix\n");
				SparsityPattern sparsity_pattern_knowns, sparsity_pattern_unknowns;
				sparsity_pattern_knowns.init(num_interior_nodes);
				sparsity_pattern_unknowns.init(num_interior_nodes);

				scalar *mem_loc;
				int ni_index, nj_index;
				for(int el = 0; el < num_elements; el++) {
					for(int ni = 0; ni < 10; ni++) {
						for(int nj = 0; nj < 10; nj++) {

							ni_index = elem[el].n[ni]->index;
							nj_index = elem[el].n[nj]->index;

							mem_loc = elem[el].get_K_alpha_element_mem_loc(ni, nj);

							/* We don't care about rows of the matrix before row num_surface_nodes */
							if(ni_index >= num_surface_nodes) {

								/* if the column is in the surface node area, add these contributions to the 'known' matrix */
								if(nj_index < num_surface_nodes) {
									sparsity_pattern_knowns.register_contribution(
															ni_index - num_surface_nodes,
															nj_index,
															mem_loc);
								}

								/* otherwise add the contribution to the 'unknowns' matrix */
								else {
									sparsity_pattern_unknowns.register_contribution(
															ni_index - num_surface_nodes,
															nj_index - num_surface_nodes,
															mem_loc);
								}
							}
						}
					}
				}

				// Use the sparsity patterns to create fixed pattern sparse matrices for use with the poisson solver
				poisson_surface_matrix = sparsity_pattern_knowns.create_sparse_matrix();
				poisson_interior_matrix = sparsity_pattern_unknowns.create_sparse_matrix();

				// Create a conjugate gradient solver for use with the 'unknowns' (interior) poisson matrix
				printf("\t\tCreating and initialising Poisson Solver...");
				poisson_solver = new CG_solver();
				poisson_solver->init(num_interior_nodes, params->epsilon2, params->max_iterations_cg);

				// Create the vector containing all values of the potential at each interior node
				phi_Omega = new scalar[num_interior_nodes];

				if(phi_Omega == NULL) {
					FFEA_ERROR_MESSG("Could not allocate memory (for phi_Omega array).\n");
				}

				for(int i = 0; i < num_interior_nodes; i++) {
					phi_Omega[i] = 0;
				}

				// Create the vector containing all values of the potential on each surface node
				phi_Gamma = new scalar[num_surface_nodes];

				if(phi_Gamma == NULL) {
					FFEA_ERROR_MESSG("Could not allocate memory (for phi_Gamma array).\n");
				}

				for(int i = 0; i < num_surface_nodes; i++) {
					phi_Gamma[i] = 0;
				}


		/*
				// Create the vector containing the smoothed out charge distribution across the elements
				q = new scalar[num_nodes];

				if(q == NULL) {
					FFEA_error_text();
					printf("Could not allocate memory (for q array).\n");
					return FFEA_ERROR;
				}

				for(int i = 0; i < num_nodes; i++) {
					q[i] = 0;
				}

				// Create the vector containing all charges on each node of the Blob
				nodal_q = new scalar[num_nodes];

				if(nodal_q == NULL) {
					FFEA_error_text();
					printf("Could not allocate memory (for nodal_q array).\n");
					return FFEA_ERROR;
				}

				for(int i = 0; i < num_nodes; i++) {
					nodal_q[i] = 0;
				}

				// Set the charges on the Blob
				for(int i = 0; i < num_nodes; i++) {
		//			nodal_q[i] = (node[i].pos.z - (2.9e-8/2.0)) * 1e35; //1e26;
					nodal_q[i] = 1e26;
		//			printf("%e %e\n", node[i].pos.z, nodal_q[i]);
				}

				// Apply mass matrix to obtain smooth extrapolated charge density across elements
				M->apply(nodal_q, q);

				for(int i = 0; i < num_nodes; i++) {
					node[i].rho = q[i];
				}
		*/

				// Create the vector containing the charge distribution across the elements
				q = new scalar[num_nodes];

				if(q == NULL) {
					FFEA_ERROR_MESSG("Could not allocate memory (for q array).\n");
				}

				for(int i = 0; i < num_nodes; i++) {
					q[i] = 0;
				}

				const scalar charge_density = 6.0e25;
				for(int n = 0; n < num_elements; n++) {
					for(int i = 0; i < 10; i++) {
						q[elem[n].n[i]->index] += charge_density * elem[n].vol_0;
					}
				}

				for(int i = 0; i < num_nodes; i++) {
					node[i].rho = q[i];
				}

		//		for(int i = 0; i < num_nodes; i++) {
		//			printf("q[%d] = %e\n", i, q[i]);
		//		}

				// Create the Right Hand Side vector for the Poisson solver
				poisson_rhs = new scalar[num_interior_nodes];

				if(poisson_rhs == NULL) {
					FFEA_ERROR_MESSG("Could not allocate memory (for poisson_rhs array).\n");
				}

				for(int i = 0; i < num_interior_nodes; i++) {
					poisson_rhs[i] = 0;
				}

		//		build_poisson_matrices_and_setup_for_solve();
		//		poisson_solver->solve(poisson_interior_matrix, phi_Omega, poisson_rhs);

				printf("\t\tdone.\n");
			}

			printf("Finished Blob initialisation.\n");

			// Return FFEA_OK to indicate "success"
			return FFEA_OK;
		}


		/*
		 * Solves the EOM on the finite element mesh, updating the node positions and velocities by one time step
		 */
		int update()
		{
			if(blob_state != FFEA_BLOB_IS_DYNAMIC) {
				return FFEA_OK;
			}

			/* some "work" variables */
			matrix3 J;			// Holds the Jacobian calculated for the *current* element being processed
			matrix3 stress;			// Holds the current stress tensor (elastic stress, with thermal fluctuations)
			vector12 du;			// Holds the force change for the current element
			int n;				// Iterates through all the elements
			int tid;			// Holds the current thread id (in parallel regions)
			int num_inversions = 0;		// Counts the number of elements that have inverted (if > 0 then simulation has failed)
	
			// Element loop
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel default(none) private(J, stress, du, tid, n) reduction(+:num_inversions)
				{
			#endif
				tid = omp_get_thread_num();

			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp for schedule(guided)
			#endif
				for(n = 0; n < num_elements; n++) {

					// calculate jacobian for this element
					elem[n].calculate_jacobian(J);

					// get the 12 derivatives of the shape functions (by inverting the jacobian)
					// and also get the element volume. The function returns an error in the
					// case of an element inverting itself (determinant changing sign since last step)
					if(elem[n].calc_shape_function_derivatives_and_volume(J) == FFEA_ERROR) {
						FFEA_error_text();
						printf("Element %d has inverted\n", n);
						num_inversions++;
					}

					// create viscosity matrix
					elem[n].create_viscosity_matrix();

					// Now build the stress tensor from the shear elastic, bulk elastic and fluctuating stress contributions
					mat3_set_zero(stress);
					elem[n].add_shear_elastic_stress(J, stress);
					elem[n].add_bulk_elastic_stress(stress);

					if(params->calc_noise == 1) {
						elem[n].add_fluctuating_stress(params, rng, stress, tid);
					}
					
					elem[n].internal_stress_mag = sqrt(mat3_double_contraction_symmetric(stress));

					// Calculate internal forces of current element (or don't, depending on solver)
					if(linear_solver != FFEA_NOMASS_CG_SOLVER) { 
						elem[n].get_element_velocity_vector(du);
						mat12_apply(elem[n].viscosity_matrix, du);
					} else {
						vec12_set_zero(du);
					}

					elem[n].apply_stress_tensor(stress, du);

					// Store the contributions to the force on each of this element's nodes (Store them on
					// the element - they will be aggregated on the actual nodes outside of this parallel region)
					elem[n].add_element_force_vector(du);

					if(params->calc_es == 1) {
						elem[n].calculate_electrostatic_forces();
					}
				}

			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				}
			#endif

			// Check if any elements have inverted
			if(num_inversions != 0) {
				if(num_inversions == 1) {
					FFEA_ERROR_MESSG("1 element has inverted since the last step. Aborting simulation.\n");
				}
				else {
					FFEA_ERROR_MESSG("%d elements have inverted since the last step. Aborting simulation.\n", num_inversions);
				}
			}
	
			// Aggregate forces on nodes from all elements
			if(aggregate_forces_and_solve() == FFEA_ERROR) {
				FFEA_ERROR_MESSG("There was a problem in 'aggregate_forces_and_solve' function.\n");
			}

			// Update node velocities and positions
			euler_integrate();

			// Linearise the 2nd order elements
			for(n = 0; n < num_elements; n++) {
				elem[n].linearise_element();
			}

			return FFEA_OK;
		}

		/*
		 * Calculates the centroid of this Blob, then translates all nodes in the Blob
		 * so that the new centroid position is at the given (x,y,z) position
		 */
		void position(scalar x, scalar y, scalar z)
		{
			int i;
			scalar centroid_x = 0, centroid_y = 0, centroid_z = 0;
			scalar dx, dy, dz;

			// Calculate centroid of [SURFACE] Blob mesh
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) private(i) reduction(+:centroid_x,centroid_y,centroid_z)
			#endif
			for(i = 0; i < num_surface_nodes; i++) {
				centroid_x += node[i].pos.x;
				centroid_y += node[i].pos.y;
				centroid_z += node[i].pos.z;
			}

			centroid_x *= (1.0/num_surface_nodes);
			centroid_y *= (1.0/num_surface_nodes);
			centroid_z *= (1.0/num_surface_nodes);

			// Calculate displacement vector required to move centroid to requested position
			dx = x - centroid_x;
			dy = y - centroid_y;
			dz = z - centroid_z;

			// Move all nodes in mesh by displacement vector
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) private(i) shared(dx, dy, dz)
			#endif
			for(i = 0; i < num_nodes; i++) {
				node[i].pos.x += dx;
				node[i].pos.y += dy;
				node[i].pos.z += dz;
			}
		}


		/*
		 * Translate the Blob by the given vector
		 */
		void move(scalar dx, scalar dy, scalar dz)
		{
			for(int i = 0; i < num_nodes; i++) {
				node[i].pos.x += dx;
				node[i].pos.y += dy;
				node[i].pos.z += dz;
			}
		}

		/*
		 * Calculates and returns the centre of mass of this Blob
		 */
		void get_CoM(vector3 *com)
		{
			com->x = 0;
			com->y = 0;
			com->z = 0;

			
			for(int n = 0; n < num_elements; n++) {
				com->x += (elem[n].mass * elem[n].n[0]->pos.x + elem[n].n[1]->pos.x + elem[n].n[2]->pos.x + elem[n].n[3]->pos.x)/4;
				com->y += (elem[n].mass * elem[n].n[0]->pos.y + elem[n].n[1]->pos.y + elem[n].n[2]->pos.y + elem[n].n[3]->pos.y)/4;
				com->z += (elem[n].mass * elem[n].n[0]->pos.z + elem[n].n[1]->pos.z + elem[n].n[2]->pos.z + elem[n].n[3]->pos.z)/4;
			}
			if(num_elements == 0) {
				return;
			} else {
				com->x /= num_elements;
				com->y /= num_elements;
				com->z /= num_elements;
			}
		}

		/*
		 * Calculates and returns the centroid of this Blob
		 */
		void get_centroid(vector3 *com)
		{
			com->x = 0;
			com->y = 0;
			com->z = 0;
			for(int n = 0; n < num_nodes; n++) {
				com->x += node[n].pos.x;
				com->y += node[n].pos.y;
				com->z += node[n].pos.z;
			}
			com->x /= num_nodes;
			com->y /= num_nodes;
			com->z /= num_nodes;
		}

		void set_rmsd_pos_0()
		{
			for(int i = 0; i < num_nodes; i++) {
				node[i].pos_0.x = node[i].pos.x;
				node[i].pos_0.y = node[i].pos.y;
				node[i].pos_0.z = node[i].pos.z;
			}
		}

		/*
		 * Writes a new node file for the case of an initially translated STATIC blob, so viewer doesn't read straight from node file
		 */
		int create_viewer_node_file(const char *node_filename, scalar scale)
		{

			FILE *out = NULL;
			char new_node_filename[] = "VIEWERNODE_";
			int i;

			// Name new node file
			strcat(new_node_filename, node_filename);
	

			//open the new node file
			if((out = fopen(new_node_filename, "w")) == NULL) {
				FFEA_FILE_ERROR_MESSG(new_node_filename);
			}
			printf("\t\tWriting to viewer nodes file: %s\n", new_node_filename);

			fprintf(out, "ffea viewer node file\n");
			fprintf(out, "num_nodes %d\n", num_nodes);
			fprintf(out, "num_surface_nodes %d\n", num_surface_nodes);
			fprintf(out, "num_interior_nodes %d\n", num_interior_nodes);

			// Write all the nodes to file
			fprintf(out, "surface nodes:\n");
			for(i = 0; i < num_surface_nodes; i++) {
				fprintf(out, "%le %le %le\n", node[i].pos.x/scale, node[i].pos.y/scale, node[i].pos.z/scale);
			}

			fprintf(out, "interior nodes:\n");
			for(i = num_surface_nodes; i < num_nodes; i++) {
				fprintf(out, "%le %le %le\n", node[i].pos.x/scale, node[i].pos.y/scale, node[i].pos.z/scale);
			}

			fclose(out);
			printf("\t\t\tWrote %d nodes from %s\n", i, new_node_filename);

			return FFEA_OK;
		}

		/*
		 * Dumps all the node positions (in order) from the node array to the given file stream.
		 */
		void write_nodes_to_file(FILE *trajectory_out)
		{
			// If this is a static blob, then don't bother printing out all the node positions (since there will be no change)
			if(blob_state == FFEA_BLOB_IS_STATIC) {
				fprintf(trajectory_out, "STATIC\n");
				return;
			} else if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
				fprintf(trajectory_out, "DYNAMIC\n");
			} else if(blob_state == FFEA_BLOB_IS_FROZEN) {
				fprintf(trajectory_out, "FROZEN\n");
			}

			fprintf(trajectory_out, "%d\n", num_nodes);
			for(int i = 0; i < num_nodes; i++) {
				fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n", node[i].pos.x, node[i].pos.y, node[i].pos.z, node[i].vel.x, node[i].vel.y, node[i].vel.z, node[i].phi, force[i].x, force[i].y, force[i].z);
			}
		}

		/*
		 * Reads the node positions from the given trajectory file stream.
		 * This is useful for restarting simulations from trajectory files.
		 */
		int read_nodes_from_file(FILE *trajectory_out)
		{
			char state_str[20];
			char *result = NULL;

			// If blob is static, don't read any nodes. Simply read the word "STATIC"
			if(blob_state == FFEA_BLOB_IS_STATIC) {
				result = fgets(state_str, 20, trajectory_out);
				if(result == NULL) {
					FFEA_ERROR_MESSG("Problem when reading 'STATIC' (expected) line in trajectory file\n")
				}
				if(strcmp(state_str, "STATIC\n") != 0) {
					FFEA_ERROR_MESSG("When restarting from trajectory file, expected to read 'STATIC', but instead found '%s...'\n", state_str)
				}
				return FFEA_OK;
			} else {
				result = fgets(state_str, 20, trajectory_out);
				if(result == NULL) {
					FFEA_ERROR_MESSG("Problem when reading state line in trajectory file\n")
				}
			}

			int read_num_nodes = 0;
			if(fscanf(trajectory_out, "%d\n", &read_num_nodes) != 1) {
				FFEA_ERROR_MESSG("Could not read number of nodes.\n")
			}
			if(read_num_nodes != num_nodes) {
				FFEA_ERROR_MESSG("(When restarting) Number of nodes in trajectory file does not match number of nodes in node file\n")
			}
			for(int i = 0; i < num_nodes; i++) {
				if(fscanf(trajectory_out, "%le %le %le %le %le %le %le %le %le %le\n", &node[i].pos.x, &node[i].pos.y, &node[i].pos.z, &node[i].vel.x, &node[i].vel.y, &node[i].vel.z, &node[i].phi, &force[i].x, &force[i].y, &force[i].z) !=10) {
					FFEA_ERROR_MESSG("(When restarting) Error reading from trajectory file, for node %d\n", i)
				}
			}
			return FFEA_OK;
		}

		/*
		 * Takes measurements of system properties: KE, PE, Centre of Mass and Angular momentum, outputing
		 * the results to the measurement file.
		 */
		void make_measurements(FILE *measurement_out, int step, vector3 *system_CoM)
		{
			// Only calculate and write out measurements if there is an open measurement output file
			if(measurement_out != NULL) {
				int n, i, j;
				scalar ke, pe, com_x, com_y, com_z, L_x, L_y, L_z;
				scalar rmsd;
				scalar temp1, temp2, temp3;
				scalar r[4][3];
				vector12 vec;

				ke = 0;	pe = 0;
				com_x = 0; com_y = 0; com_z = 0;
				L_x = 0; L_y = 0; L_z = 0;
				#ifdef FFEA_PARALLEL_WITHIN_BLOB
					#pragma omp parallel for default(none) reduction(+:ke, pe, com_x, com_y, com_z) private(n, vec, temp1, temp2)
				#endif
				for(n = 0; n < num_elements; n++) {

					if(linear_solver != FFEA_NOMASS_CG_SOLVER) {

						/*
			 			 * Kinetic energy contribution:
			 			 */
						
						// Read the u vector for this element
						elem[n].get_element_velocity_vector(vec);

						// Apply the mass matrix
						elem[n].apply_element_mass_matrix(vec);

						// Dot u with M.u to get the contribution to the kinetic energy
						ke +=   elem[n].n[0]->vel.x * vec[0] +
							elem[n].n[1]->vel.x * vec[1] +
							elem[n].n[2]->vel.x * vec[2] +
							elem[n].n[3]->vel.x * vec[3] +
							elem[n].n[0]->vel.y * vec[4] +
							elem[n].n[1]->vel.y * vec[5] +
							elem[n].n[2]->vel.y * vec[6] +
							elem[n].n[3]->vel.y * vec[7] +
							elem[n].n[0]->vel.z * vec[8] +
							elem[n].n[1]->vel.z * vec[9] +
							elem[n].n[2]->vel.z * vec[10] +
							elem[n].n[3]->vel.z * vec[11];

						/*
			 			 * Centre of Mass contribution:
			 			 */
						com_x += elem[n].mass * (elem[n].n[0]->pos.x + elem[n].n[1]->pos.x + elem[n].n[2]->pos.x + elem[n].n[3]->pos.x)/4;
						com_y += elem[n].mass * (elem[n].n[0]->pos.y + elem[n].n[1]->pos.y + elem[n].n[2]->pos.y + elem[n].n[3]->pos.y)/4;
						com_z += elem[n].mass * (elem[n].n[0]->pos.z + elem[n].n[1]->pos.z + elem[n].n[2]->pos.z + elem[n].n[3]->pos.z)/4;
					}
	
					/*
			 		 * Potential energy contribution:
			 		 */

					// Old
					/*temp2 = elem[n].G/elem[n].E;
					temp1 = (elem[n].vol/elem[n].vol_0) - (1 + temp2);
					pe += elem[n].vol_0 * (
				  		elem[n].G * ( mat3_double_contraction_symmetric(elem[n].F_ij) - 3 )
						+ elem[n].E * ( temp1 * temp1 - temp2 * temp2 ));*/

					// New 
					scalar C = elem[n].E - elem[n].G * 2.0/3.0;
					temp1 = elem[n].vol/elem[n].vol_0;
					pe += elem[n].vol_0 * (
						elem[n].G * ( mat3_double_contraction(elem[n].F_ij) - 3 )
						+ 0.5 * C * (pow(temp1, 2) - 1)
						- ((2 * elem[n].G) + C) * log(temp1)); 

				}

				// And don't forget to multiply by a half
				ke *= .5;
				pe *= .5;

				if(linear_solver != FFEA_NOMASS_CG_SOLVER) {

					// Get the centre of mass from the calculated totals
					if(mass <= 0.0) {
						com_x = 0.0;
						com_y = 0.0;
						com_z = 0.0;
					} else {
						com_x *= 1.0/mass;
						com_y *= 1.0/mass;
						com_z *= 1.0/mass;
					}	
				
	
					/* Calculate angular momentum */
					#ifdef FFEA_PARALLEL_WITHIN_BLOB
						#if defined(__GNUC__)
							#pragma omp parallel for default(none) reduction(+:L_x,L_y,L_z) shared(system_CoM) private(n, r, temp1, temp2, temp3, i, j)
						#elif defined(__ICC)
							#pragma omp parallel for default(none) reduction(+:L_x,L_y,L_z) shared(system_CoM) private(n, r, temp1, temp2, temp3, i, j)
						#endif
					#endif
					for(n = 0; n < num_elements; n++) {
	
						// mass matrix
						const matrix4 M =      {{.1, .05, .05, .05},
									{.05, .1, .05, .05},
									{.05, .05, .1, .05},
									{.05, .05, .05, .1}};
	
						// Find the separation vectors for this element
						for(i = 0; i < 4; i++) {
							r[i][0] = elem[n].n[i]->pos.x - system_CoM->x;
							r[i][1] = elem[n].n[i]->pos.y - system_CoM->y;
							r[i][2] = elem[n].n[i]->pos.z - system_CoM->z;
						}

						// Calculate contribution to angular momentum from this element
						temp1 = 0; temp2 = 0; temp3 = 0;
						for(i = 0; i < 4; i++) {
							for(j = 0; j < 4; j++) {
								temp1 += M[i][j] * (r[j][1] * elem[n].n[i]->vel.z - r[j][2] * elem[n].n[i]->vel.y);
								temp2 += M[i][j] * (r[j][2] * elem[n].n[i]->vel.x - r[j][0] * elem[n].n[i]->vel.z);
								temp3 += M[i][j] * (r[j][0] * elem[n].n[i]->vel.y - r[j][1] * elem[n].n[i]->vel.x);
							}
						}
	
						// Add this contribution to the sum
						L_x += elem[n].mass * temp1;
						L_y += elem[n].mass * temp2;
						L_z += elem[n].mass * temp3;
					}
				}

				/* Calculate RMSD value for this configuration */
				rmsd = 0;
				#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) reduction(+:rmsd) private(i, temp1, temp2, temp3)
				#endif
				for(i = 0; i < num_nodes; i++) {
					temp1 = node[i].pos.x - node[i].pos_0.x;
					temp2 = node[i].pos.y - node[i].pos_0.y;
					temp3 = node[i].pos.z - node[i].pos_0.z;
					rmsd += temp1*temp1 + temp2*temp2 + temp3*temp3;
				}
				rmsd = sqrt(rmsd/num_nodes);

				/* Count cumulative area of faces engaging in VdW_xz on this blob */
				/* For those that contribute, calculate the VdW energy contribution */
				scalar vdw_xz_area = 0;
				scalar vdw_xz_energy = 0;
				vector3 tot_surf_force = {0,0,0};
				scalar vdw_xz_force_mag = 0.0;
				scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
				lj_matrix->get_LJ_params(0, 0, &vdw_eps, &vdw_r_eq);

				for(i = 0; i < num_surface_faces; i++) {
					if(surface[i].vdw_xz_interaction_flag == true) {
						vdw_xz_area += surface[i].area;
						vdw_xz_energy += surface[i].get_vdw_energy_with_xz_plane(params->es_h * (1.0/params->kappa) * params->es_N_y, vdw_r_eq, vdw_eps);
					}
				}

				/* Calculate the total VdW_xz force on the Blob */
				for(n = 0; n < num_surface_faces; n++) {
					if(surface[i].vdw_xz_interaction_flag == true) {
						for(int i = 0; i < 3; i++) {
							tot_surf_force.x += surface[n].force[i].x;
							tot_surf_force.y += surface[n].force[i].y;
							tot_surf_force.z += surface[n].force[i].z;
						}
					}
				}
				vdw_xz_force_mag = mag(&tot_surf_force);

				// print out all measurements to file
				fprintf(measurement_out, "%d\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\n", step, ke, pe, com_x, com_y, com_z, L_x, L_y, L_z, rmsd, vdw_xz_area, vdw_xz_force_mag, vdw_xz_energy);
			}
		}

		void make_stress_measurements(FILE *stress_out, int blob_number)
		{
			int n;
			if(stress_out != NULL) {
				fprintf(stress_out, "blob\t%d\n", blob_number);
				for(n = 0; n < num_elements; n++) {
					fprintf(stress_out, "%e\n", elem[n].internal_stress_mag);
				}
				fprintf(stress_out, "\n");
			}
		}

		/*
		 * Get the centroid of all faces on the blob surface
		 */
		void calc_centroids_and_normals_of_all_faces()
		{
			int i;
			for(i = 0; i < num_surface_faces; i++)
				surface[i].calc_area_normal_centroid();
		}

		/*
		 *
		 */
		int get_num_faces()
		{
			return num_surface_faces;
		}

		/*
		 * Return pointer to the ith Face of this Blob's surface
		 */
		Face *get_face(int i)
		{
			if(surface[i].is_vdw_active() == true) {
				return &surface[i];
			} else {
				return NULL;
			}
		}

		scalar get_vdw_area()
		{
			scalar total_vdw_area = 0.0;
			for(int i = 0; i < get_num_faces(); ++i) {
				if(surface[i].is_vdw_active() == true) {
					total_vdw_area += surface[i].area;
				}
			}
			return total_vdw_area;
		}
//		/*
//		 *
//		 */
//		int get_num_surface_elements();


		/*
		 * Solves the poisson equation inside this Blob for a given fixed surface potential, and calculates the surface flux out of the protein
		 */
		int solve_poisson(scalar *phi_gamma_IN, scalar *J_Gamma_OUT)
		{
			if(num_interior_nodes > 0) {
				/* Convert the given potential on the surface faces into the potential on each node */
				for(int i = 0; i < num_surface_nodes; i++) {
					phi_Gamma[i] = 0;
				}

				for(int i = 0; i < num_surface_faces; i++) {
					for(int j = 0; j < 3; j++) {
						phi_Gamma[surface[i].n[j]->index] += phi_gamma_IN[i];
					}
				}

				for(int i = 0; i < num_surface_nodes; i++) {
					phi_Gamma[i] /= num_contributing_faces[i];
		//			printf("phi_Gamma[%d] = %e\n", i, phi_Gamma[i]);
				}

				/* Calculate the RHS vector */
				poisson_surface_matrix->apply(phi_Gamma, poisson_rhs);
				for(int i = 0; i < num_interior_nodes; i++) {
					poisson_rhs[i] = q[i + num_surface_nodes] - poisson_rhs[i];
		//			printf("poisson_rhs[%d] = %e\n", i, poisson_rhs[i]);
				}

				poisson_solver->solve(poisson_interior_matrix, phi_Omega, poisson_rhs);

		//		printf("int:\n");
		//		poisson_interior_matrix->print_dense();

				for(int n = 0; n < num_surface_nodes; n++)
					node[n].phi = phi_Gamma[n];

				for(int n = 0; n < num_interior_nodes; n++) {
					node[n + num_surface_nodes].phi = phi_Omega[n];
		//			printf("phi_Omega[%d] = %e\n", n, phi_Omega[n]);
				}

				for(int i = 0; i < num_surface_faces; i++) {
					J_Gamma_OUT[i] = surface[i].get_normal_flux();
				}
			}

			return FFEA_OK;
		}


		/*
		 * Set all forces on the blob to zero
		 */
		void zero_force()
		{
			for(int i = 0; i < num_elements; i++) {
				elem[i].zero_force();
			}
			for(int i = 0; i < num_surface_faces; i++) {
				surface[i].zero_force();
			}
		}


		/*
		 * Set all nodes on the Blob to the given velocity vector
		 */
		void velocity_all(scalar vel_x, scalar vel_y, scalar vel_z)
		{
			int i;
			for(i = 0; i < num_nodes; i++) {
				node[i].vel.x = vel_x;
				node[i].vel.y = vel_y;
				node[i].vel.z = vel_z;
			}
		}

		/*
		 * Constructs the Poisson matrix for this Blob.
		 */
		void build_poisson_matrices()
		{
			if(num_interior_nodes > 0) {
		/*
				scalar K[num_nodes*num_nodes];
				for(int i = 0; i < num_nodes * num_nodes; i++) {
					K[i] = 0;
				}
		*/

				/* Calculate the K_alpha matrices for each element (the diffusion matrix * epsilon * volume) */
				for(int n = 0; n < num_elements; n++) {
					elem[n].calculate_K_alpha();
		//			elem[n].add_K_alpha(K, num_nodes);
				}


		/*
				printf("BLOB Full K:\n");
				int c = 0;
				for(int i = 0; i < num_nodes; i++) {
					for(int j = 0; j < num_nodes; j++) {
						printf("%+e ", K[c]);
						c++;
					}
					printf("\n");
				}
		*/
				/* Construct the poisson matrices for the current blob (based on diffusion matrices of elements) */
				poisson_surface_matrix->build();
				poisson_interior_matrix->build();
			}
		}

//		int elements_are_connected(int e1, int e2);

		/*
		 * Returns the total mass of this Blob.
		 */
		scalar get_mass()
		{
			return mass;
		}

		/*
		 * Applies the WALL_TYPE_HARD boundary conditions simply by finding all nodes that have "passed through" the wall,
		 * then zeroing any component of those nodes' velocities that points into the wall. This allows the Blob to "wobble"
		 * its way back out of the wall (and effectively prevents any penetration larger than dt * largest velocity,
		 * generally very small for sensible dt).
		 */
		void enforce_box_boundaries(vector3 *box_dim)
		{
			if(params->wall_x_1 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.x < 0) {
						if(node[i].vel.x < 0) {
							node[i].vel.x = 0;
						}
					}
				}
			}
			if(params->wall_x_2 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.x > box_dim->x) {
						if(node[i].vel.x > 0) {
							node[i].vel.x = 0;
						}
					}
				}
			}
			if(params->wall_y_1 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.y < 0) {
						if(node[i].vel.y < 0) {
							node[i].vel.y = 0;
						}
					}
				}
			}
			if(params->wall_y_2 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.y > box_dim->y) {
						if(node[i].vel.y > 0) {
							node[i].vel.y = 0;
						}
					}
				}
			}
			if(params->wall_z_1 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.z < 0) {
						if(node[i].vel.z < 0) {
							node[i].vel.z = 0;
						}
					}
				}
			}
			if(params->wall_z_2 == WALL_TYPE_HARD) {
				for(int i = 0; i < num_surface_nodes; i++) {
					if(node[i].pos.z > box_dim->z) {
						if(node[i].vel.z > 0) {
							node[i].vel.z = 0;
						}
					}
				}
			}
		}

		/*
		 * Set the interaction flag for all faces of this Blob back to false (i.e. not interacting)
		 */
		void reset_all_faces()
		{
			for(int i = 0; i < num_surface_faces; i++) {
				surface[i].set_vdw_xz_interaction_flag(false);
			}
		}


	private:

		/* Total number of nodes in Blob */
		int num_nodes;

		/* Total number of elements in Blob */
		int num_elements;

		/* Total number of surface elements in Blob */
		int num_surface_elements;

		/* Total number of interior elements in Blob */
		int num_interior_elements;

		/* Total number of faces on Blob surface */
		int num_surface_faces;

		/* Total number of nodes on Blob surface */
		int num_surface_nodes;

		/* Total number of nodes on Blob interior */
		int num_interior_nodes;

		/* Number of 'pinned' nodes (nodes which are not able to move, removing degrees of freedom from system) */
		int num_pinned_nodes;

		/* Whether this Blob is DYNAMIC (movable; dynamics simulated) or STATIC (fixed; no simulation of dynamics; Blob is a perfectly solid object fixed in space)*/
		int blob_state;

		/* Total mass of Blob */
		scalar mass;

		/* Total vdw energy between blobs */
		scalar vdw_bb_energy;

		/* Array of nodes */
		mesh_node *node;

		/* Array of elements */
		tetra_element_linear *elem;

		/* Array of surface faces */
		Face *surface;

		/* List of fixed ('pinned') nodes */
		int *pinned_nodes_list;

		/* A pointer to a class containing simulation parameters, such as the time step, dt */
		SimulationParams *params;

		/* pointer to the vdw forcefield parameters (for some energy calcs) */
		LJ_matrix *lj_matrix;

		/* A pointer to whatever Solver is being used for this Blob (eg SparseSubstitutionSolver
		 * or ConjugateGradientSolver). The Solver solves the equation Mx = f where M is the
		 * mass matrix of this Blob and f is the force vector.
		 */
		Solver *solver;

		/* Remember what type of solver we are using */
		int linear_solver;

		/* The Blob force vector (an array of the force on every node) */
		vector3 *force;

		/* The array of random number generators (needed for parallel runs) */
		MTRand *rng;

		CG_solver *poisson_solver;
		SparseMatrixFixedPattern *poisson_surface_matrix;
		SparseMatrixFixedPattern *poisson_interior_matrix;
		scalar *phi_Omega;
		scalar *phi_Gamma;
		scalar *q;
		scalar *nodal_q;
		scalar *poisson_rhs;

		int *num_contributing_faces;

		connectivity_entry *element_connectivity_table;

		/*
		 *
		 */
		SparseMatrixFixedPattern *M;


		/*
 		* Opens and reads the given 'ffea node file', extracting all the nodes for this Blob.
 		* Records how many of these are surface nodes and how many are interior nodes.
 		*/
		int load_nodes(const char *node_filename, scalar scale)
		{
			FILE *in = NULL;
			int i;
			double x, y, z;
			const int max_line_size = 50;
			char line[max_line_size];

			// open the node file
			if((in = fopen(node_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(node_filename)
			}
			printf("\t\tReading in nodes file: %s\n", node_filename);

			// first line should be the file type "ffea node file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of node file\n")
			}
			if(strcmp(line, "walrus node file\n") != 0 && strcmp(line, "ffea node file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea node file' (read '%s') \n", line)
			}

			// read in the number of nodes in the file
			if(fscanf(in, "num_nodes %d\n", &num_nodes) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of nodes\n")
			}
			printf("\t\t\tNumber of nodes = %d\n", num_nodes);

			// read in the number of surface nodes in the file
			if(fscanf(in, "num_surface_nodes %d\n", &num_surface_nodes) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of surface nodes\n")
			}
			printf("\t\t\tNumber of surface nodes = %d\n", num_surface_nodes);

			// read in the number of interior nodes in the file
			if(fscanf(in, "num_interior_nodes %d\n", &num_interior_nodes) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of interior nodes\n")
			}
			printf("\t\t\tNumber of interior nodes = %d\n", num_interior_nodes);

			// Allocate the memory for all these nodes
			node = new mesh_node[num_nodes];
			if(node == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Unable to allocate memory for nodes array.\n")
			}

			// Check for "surface nodes:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'surface nodes:' line\n")
			}
			if(strcmp(line, "surface nodes:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'surface nodes:' line (found '%s' instead)\n", line)
			}

			// Read in all the surface nodes from file
			for(i = 0; i < num_nodes; i++) {

				if(i == num_surface_nodes) {
					// Check for "interior nodes:" line
					if(fgets (line, max_line_size, in) == NULL) {
						fclose(in);
						FFEA_ERROR_MESSG("Error when looking for 'interior nodes:' line\n")
					}
					if(strcmp(line, "interior nodes:\n") != 0) {
						fclose(in);
						FFEA_ERROR_MESSG("Could not find 'interior nodes:' line (found '%s' instead)\n", line)
					}
				}

				if(fscanf(in, "%le %le %le\n", &x, &y, &z) != 3) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from nodes file at node %d\n", i)
				} else {
					node[i].pos.x = scale * (scalar) x;
					node[i].pos.y = scale * (scalar) y;
					node[i].pos.z = scale * (scalar) z;

					node[i].vel.x = 0;
					node[i].vel.y = 0;
					node[i].vel.z = 0;

					node[i].index = i;
				}
			}

			fclose(in);
			printf("\t\t\tRead %d nodes from %s\n", i, node_filename);

			return FFEA_OK;
		}

		/*
 		* Opens and reads the given 'ffea topology file', extracting all the elements for this Blob.
 		* Records how many of these are surface elements and how many are interior elements.
 		*/
		int load_topology(const char *topology_filename)
		{
			FILE *in = NULL;
			int i, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
			const int max_line_size = 50;
			char line[max_line_size];

			// Now open the topology file
			if((in = fopen(topology_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(topology_filename)
			}
			printf("\t\tReading in topology file: %s\n", topology_filename);

			// first line should be the file type "ffea topology file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of topology file\n")
			}
			if(strcmp(line, "walrus topology file\n") != 0 && strcmp(line, "ffea topology file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea topology file' (read '%s') \n", line)
			}

			// read in the total number of elements in the file
			if(fscanf(in, "num_elements %d\n", &num_elements) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of elements\n")
			}
			printf("\t\t\tNumber of elements = %d\n", num_elements);

			// read in the number of surface elements in the file
			if(fscanf(in, "num_surface_elements %d\n", &num_surface_elements) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of surface elements\n")
			}
			printf("\t\t\tNumber of surface elements = %d\n", num_surface_elements);

			// read in the number of interior elements in the file
			if(fscanf(in, "num_interior_elements %d\n", &num_interior_elements) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of interior elements\n")
			}
			printf("\t\t\tNumber of interior elements = %d\n", num_interior_elements);

			// Allocate the memory for all these elements
			elem = new tetra_element_linear[num_elements];
			if(elem == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Unable to allocate memory for element array.\n")
			}

			// Check for "surface elements:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'surface elements:' line\n")
			}
			if(strcmp(line, "surface elements:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'surface elements:' line (found '%s' instead)\n", line)
			}

			// Read in all the elements from file
			for(i = 0; i < num_elements; i++) {
				if(i == num_surface_elements) {
					// Check for "interior elements:" line
					if(fgets (line, max_line_size, in) == NULL) {
						fclose(in);
						FFEA_ERROR_MESSG("Error when looking for 'interior elements:' line\n")
					}
					if(strcmp(line, "interior elements:\n") != 0) {
						fclose(in);
						FFEA_ERROR_MESSG("Could not find 'interior elements:' line (found '%s' instead)\n", line)
					}
				}

				if(fscanf(in, "%d %d %d %d %d %d %d %d %d %d\n", &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9, &n10) != 10) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from elements file at element %d\n", i)
				} else {
					// check that none of these reference nodes outside of the node array
					if(	n1 < 0 || n1 >= num_nodes ||
						n2 < 0 || n2 >= num_nodes ||
						n3 < 0 || n3 >= num_nodes ||
						n4 < 0 || n4 >= num_nodes ||
						n5 < 0 || n5 >= num_nodes ||
						n6 < 0 || n6 >= num_nodes ||
						n7 < 0 || n7 >= num_nodes ||
						n8 < 0 || n8 >= num_nodes ||
						n9 < 0 || n9 >= num_nodes ||
						n10 < 0 || n10 >= num_nodes)
					{
						fclose(in);
						FFEA_ERROR_MESSG("Error: Element %d references an out of bounds node index\n", i)
					}

					elem[i].n[0] = &node[n1];
					elem[i].n[1] = &node[n2];
					elem[i].n[2] = &node[n3];
					elem[i].n[3] = &node[n4];
					elem[i].n[4] = &node[n5];
					elem[i].n[5] = &node[n6];
					elem[i].n[6] = &node[n7];
					elem[i].n[7] = &node[n8];
					elem[i].n[8] = &node[n9];
					elem[i].n[9] = &node[n10];

					elem[i].daddy_blob = this;
					elem[i].index = i;
				}
			}

			fclose(in);

			if(i==1)
				printf("\t\t\tRead 1 element from %s\n", topology_filename);
			else
				printf("\t\t\tRead %d elements from %s\n", i, topology_filename);

			return FFEA_OK;
		}

		/*
 		* Opens and reads the given 'ffea surface file', extracting all the faces for this Blob.
 		*/
		int load_surface(const char *surface_filename)
		{
			FILE *in = NULL;
			int i, n1, n2, n3, element;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(surface_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(surface_filename)
			}
			printf("\t\tReading in surface file: %s\n", surface_filename);

			// first line should be the file type "ffea surface file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of surface file\n")
			}
			if(strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea surface file' (read '%s') \n", line)
			}

			// read in the number of faces in the file
			if(fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of faces\n")
			}
			printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

			// Allocate the memory for all these faces
			surface = new Face[num_surface_faces];

			// Check for "faces:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'faces:' line\n")
			}
			if(strcmp(line, "faces:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'faces:' line (found '%s' instead)\n", line)
			}

			// Read in all the faces from file
			scalar smallest_A = INFINITY;
			for(i = 0; i < num_surface_faces; i++) {
				if(fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from surface file at face %d. There should be 4 space separated integers. \n", i);
				} else {
					// check that none of these reference nodes outside of the node array
					if(	n1 < 0 || n1 >= num_nodes ||
						n2 < 0 || n2 >= num_nodes ||
						n3 < 0 || n3 >= num_nodes)
					{
						fclose(in);
						FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds node index\n", i);
					} else if(element < 0 || element >= num_elements) {
						fclose(in);
						FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds element index\n", i);
					}


					int n1_el = elem[element].what_node_is_this(n1), n2_el = elem[element].what_node_is_this(n2), n3_el = elem[element].what_node_is_this(n3);

					SecondOrderFunctions::stu n1_stu = {SecondOrderFunctions::stu_lookup[n1_el].s, SecondOrderFunctions::stu_lookup[n1_el].t, SecondOrderFunctions::stu_lookup[n1_el].u};
					SecondOrderFunctions::stu n2_stu = {SecondOrderFunctions::stu_lookup[n2_el].s, SecondOrderFunctions::stu_lookup[n2_el].t, SecondOrderFunctions::stu_lookup[n2_el].u};
					SecondOrderFunctions::stu n3_stu = {SecondOrderFunctions::stu_lookup[n3_el].s, SecondOrderFunctions::stu_lookup[n3_el].t, SecondOrderFunctions::stu_lookup[n3_el].u};

					SecondOrderFunctions::stu centroid_stu =
						{
							(n1_stu.s + n2_stu.s + n3_stu.s)/3.0,
							(n1_stu.t + n2_stu.t + n3_stu.t)/3.0,
							(n1_stu.u + n2_stu.u + n3_stu.u)/3.0
						};

					surface[i].init(&elem[element], &node[n1], &node[n2], &node[n3], centroid_stu, this);

					if(surface[i].area_0 < smallest_A) {
						smallest_A = surface[i].area_0;
					}
				}
			}
			fclose(in);

			printf("\t\t\tSmallest Face Area = %e\n", smallest_A);
			if(i==1)
				printf("\t\t\tRead 1 surface face from %s\n", surface_filename);
			else
				printf("\t\t\tRead %d surface faces from %s\n", i, surface_filename);

			return FFEA_OK;
		}


		/*
 		* Opens and reads the given 'ffea surface file', extracting all the faces for this Blob, ignoring topology.
 		*/
		int load_surface_no_topology(const char *surface_filename)
		{
			FILE *in = NULL;
			int i, n1, n2, n3, element;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(surface_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(surface_filename)
			}
			printf("\t\tReading in surface file: %s\n", surface_filename);

			// first line should be the file type "ffea surface file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of surface file\n")
			}
			if(strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea surface file' (read '%s') \n", line)
			}

			// read in the number of faces in the file
			if(fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of faces\n")
			}
			printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

			// Allocate the memory for all these faces
			surface = new Face[num_surface_faces];

			// Check for "faces:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'faces:' line\n")
			}
			if(strcmp(line, "faces:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'faces:' line (found '%s' instead)\n", line)
			}

			// Read in all the faces from file
			scalar smallest_A = INFINITY;
			for(i = 0; i < num_surface_faces; i++) {
				if(fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from surface file at face %d. There should be 4 space separated integers. \n", i);
				} else {
					// check that none of these reference nodes outside of the node array
					if(	n1 < 0 || n1 >= num_nodes ||
						n2 < 0 || n2 >= num_nodes ||
						n3 < 0 || n3 >= num_nodes)
					{
						fclose(in);
						FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds node index\n", i);
					}

					surface[i].init(&node[n1], &node[n2], &node[n3], this);

					if(surface[i].area_0 < smallest_A) {
						smallest_A = surface[i].area_0;
					}
				}
			}

			fclose(in);

			printf("\t\t\tSmallest Face Area = %e\n", smallest_A);
			if(i==1)
				printf("\t\t\tRead 1 surface face from %s\n", surface_filename);
			else
				printf("\t\t\tRead %d surface faces from %s\n", i, surface_filename);

			return FFEA_OK;
		}

		/*
 		 * Opens and reads the given 'ffea material params file', extracting the material parameters for each element.
 		 */
		int load_material_params(const char *material_params_filename)
		{
			FILE *in = NULL;
			int num_material_elements;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(material_params_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(material_params_filename)
			}
			printf("\t\tReading in material parameters file: %s\n", material_params_filename);

			// first line should be the file type "ffea material params file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of material params file\n")
			}
			if(strcmp(line, "walrus material params file\n") != 0 && strcmp(line, "ffea material params file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea material params file' (read '%s') \n", line)
			}

			// read in the number of elements in the file
			if(fscanf(in, "num_elements %d\n", &num_material_elements) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of elements\n")
			}
			printf("\t\t\tNumber of elements in material params file = %d\n", num_material_elements);

			// Check that we have same number of elements in material params file as in topology file
			if(num_material_elements != num_elements) {
				fclose(in);
				FFEA_ERROR_MESSG("Number of elements in material params file (%d) does not match number of elements in topology file (%d)\n", num_material_elements, num_elements)
			}

			// Set the material parameters for each element in the Blob
			scalar density=0.0, shear_visc=0.0, bulk_visc=0.0, shear_mod=0.0, bulk_mod=0.0, dielectric=0.0;
			int i;
			for(i = 0; i < num_elements; i++) {
				if(fscanf(in, "%le %le %le %le %le %le\n", &density, &shear_visc, &bulk_visc, &shear_mod, &bulk_mod, &dielectric) != 6) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from material params file at element %d. There should be 6 space separated real values (density, shear_visc, bulk_visc, shear_mod, bulk_mod, dielectric).\n", i);
				}
				elem[i].rho = density;
				elem[i].A = shear_visc;
				elem[i].B = bulk_visc - (2.0/3.0) * shear_visc; // Code uses second coefficient of viscosity
				elem[i].G = shear_mod;
				elem[i].E = bulk_mod;
				elem[i].dielectric = dielectric;
			}

			fclose(in);

			printf("\t\t\tRead %d element material params from %s\n", i, material_params_filename);

			return FFEA_OK;
		}

		/*
 		* Opens and reads the given 'ffea stokes params file', extracting the stokes radii for each node in the Blob.
 		*/
		int load_stokes_params(const char *stokes_filename, scalar scale)
		{
			FILE *in = NULL;
			int num_stokes_nodes;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(stokes_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(stokes_filename)
			}
			printf("\t\tReading in material parameters file: %s\n", stokes_filename);

			// first line should be the file type "ffea stokes radii file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of stokes radii file\n")
			}
			if(strcmp(line, "walrus stokes radii file\n") != 0 && strcmp(line, "ffea stokes radii file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea stokes radii file' (read '%s') \n", line)
			}

			// read in the number of nodes in the file
			if(fscanf(in, "num_nodes %d\n", &num_stokes_nodes) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of nodes\n")
			}
			printf("\t\t\tNumber of nodes in stokes radii file = %d\n", num_stokes_nodes);

			// Check that we have same number of nodes in stokes radii file as in nodes file
			if(num_stokes_nodes != num_nodes) {
				fclose(in);
				FFEA_ERROR_MESSG("Number of nodes in stokes radii file (%d) does not match number of nodes in nodes file (%d)\n", num_stokes_nodes, num_nodes)
			}

			// Set the stokes radius for each node in the Blob
			scalar stokes_radius=0.0;
			int i, check = 0;
			for(i = 0; i < num_nodes; i++) {
				if(fscanf(in, "%le\n", &stokes_radius) != 1) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from stokes radii file at node %d. There should be 1 real value (stokes radius).\n", i);
				}
				node[i].stokes_radius = stokes_radius * scale;
				if(node[i].stokes_radius < 1e-12 && node[i].stokes_radius > 0.0 && check == 0) {
					char finish[1];
					int done = 0;
					printf("WARNING. Stokes Radius on node %d in this Blob is very small, %e.\nStokes Radius is scaled by same factor as node positions, so specify in same units.\n", i, node[i].stokes_radius);
					printf("Would you like to continue (y or n)?:");
					while (done == 0) {
						scanf("%s", finish);
						if(strcmp(finish, "y") == 0) {
							done = 1;
							check = 1;
						} else if (strcmp(finish, "n") == 0) {
							return FFEA_ERROR;
						} else {
							printf("Please enter y or n:");
						}
						
					}
				}
			}

			fclose(in);

			printf("\t\t\tRead %d stokes radii from %s\n", i, stokes_filename);

			return FFEA_OK;
		}


		/*
 		* Opens and reads the given 'ffea vdw file', extracting all the van der waals species for each face of this Blob.
 		*/
		int load_vdw(const char *vdw_filename, int num_vdw_face_types)
		{
			FILE *in = NULL;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(vdw_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(vdw_filename)
			}
			printf("\t\tReading in Van der Waals file: %s\n", vdw_filename);

			// first line should be the file type "ffea vdw file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of VdW file\n")
			}
			if(strcmp(line, "walrus vdw file\n") != 0 && strcmp(line, "ffea vdw file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea vdw file' (read '%s') \n", line)
			}

			// read in the number of faces in the file
			int num_vdw_faces = 0;
			if(fscanf(in, "num_faces %d\n", &num_vdw_faces) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of faces\n")
			}
			printf("\t\t\tNumber of faces = %d\n", num_vdw_faces);

			if(num_vdw_faces != num_surface_faces) {
				FFEA_ERROR_MESSG("Number of faces specified in Van der Waals file (%d) does not agree with number in surface file (%d)\n", num_vdw_faces, num_surface_faces)
			}

			// Check for "vdw params:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'vdw params:' line\n")
			}
			if(strcmp(line, "vdw params:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'vdw params:' line (found '%s' instead)\n", line)
			}

			// Read in all the vdw parameters from the file, assigning them to the appropriate faces
			int i;
			int vdw_type = 0.0;
			for(i = 0; i < num_surface_faces; i++) {
				if(fscanf(in, "%d\n", &vdw_type) != 1) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from vdw file at face %d. There should be 1 integer denoting vdw face species (-1 - unreactive). \n", i);
				} else {
					if(vdw_type > num_vdw_face_types - 1) {
						FFEA_ERROR_MESSG("Error reading from vdw file at face %d. The given vdw face type (%d) is higher than that allowed by the vdw forcefield params file (%d). \n", i, vdw_type, num_vdw_face_types - 1);
					}
					surface[i].set_vdw_interaction_type(vdw_type);
				}
			}

			fclose(in);

			printf("\t\t\tRead %d vdw faces from %s\n", i, vdw_filename);

			return FFEA_OK;
		}

		/*
 		* Opens and reads the given 'ffea pinned nodes file', extracting all the faces for this Blob.
 		*/
		int load_pinned_nodes(const char *pin_filename)
		{
			FILE *in = NULL;
			int i, pn_index;
			const int max_line_size = 50;
			char line[max_line_size];

			if((in = fopen(pin_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(pin_filename)
			}
			printf("\t\tReading in pinned nodes file: %s\n", pin_filename);

			// first line should be the file type "ffea pinned nodes file"
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading first line of pin file\n")
			}
			if(strcmp(line, "walrus pinned nodes file\n") != 0 && strcmp(line, "ffea pinned nodes file\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("This is not a 'ffea pinned nodes file' (read '%s') \n", line)
			}

			// read in the number of pinned node indices in the file
			if(fscanf(in, "num_pinned_nodes %d\n", &num_pinned_nodes) != 1) {
				fclose(in);
				FFEA_ERROR_MESSG("Error reading number of pinned nodes\n")
			}
			printf("\t\t\tNumber of pinned nodes = %d\n", num_pinned_nodes);

			// Allocate the memory for the list of pinned node indices
			pinned_nodes_list = new int[num_pinned_nodes];

			// Check for "pinned nodes:" line
			if(fgets (line, max_line_size, in) == NULL) {
				fclose(in);
				FFEA_ERROR_MESSG("Error when looking for 'pinned nodes:' line\n")
			}
			if(strcmp(line, "pinned nodes:\n") != 0) {
				fclose(in);
				FFEA_ERROR_MESSG("Could not find 'pinned nodes:' line (found '%s' instead)\n", line)
			}

			// Read in all the pinned node indices from file
			for(i = 0; i < num_pinned_nodes; i++) {
				if(fscanf(in, "%d\n", &pn_index) != 1) {
					fclose(in);
					FFEA_ERROR_MESSG("Error reading from pinned node file at face %d. There should be 1 integer per line. \n", i);
				} else {
					// check that this does not reference nodes outside of the node array
					if(pn_index < 0 || pn_index >= num_nodes) {
						fclose(in);
						FFEA_ERROR_MESSG("Error: Pinned node %d references an out of bounds node index\n", i);
					}
					pinned_nodes_list[i] = pn_index;
				}
			}

			fclose(in);

			if(i==1)
				printf("\t\t\tRead 1 pinned node index from %s\n", pin_filename);
			else
				printf("\t\t\tRead %d pinned node indices from %s\n", i, pin_filename);

			return FFEA_OK;
		}
		
		
		/*
		 * Calculate some quantities such as the rest jacobian, rest volume etc.
		 */
		void calc_rest_state_info()
		{
			// Calculate some rest state information for each element:
			// Calculate the inverse jacobian (needed for gradient deformation tensor
			// calc.) and use the determinant to calculate the element's rest volume
			// (needed for volumetric spring calc.). Also, get the diffusion matrix info
			// for constructing the preliminary poisson solver matrix.
			matrix3 J;
			scalar min_vol = INFINITY, temp;
			int min_vol_elem = 0;
			mass = 0;
			scalar total_vol = 0;
			for(int i = 0; i < num_elements; i++) {
				// Get jacobian matrix for this element
				elem[i].calculate_jacobian(J);

				// Get the inverse jacobian matrix
				mat3_invert(J, elem[i].J_inv_0, &temp);

				// get the 12 derivatives of the shape functions (by inverting the jacobian)
				// and also get the element volume
				elem[i].calc_shape_function_derivatives_and_volume(J);
				elem[i].vol_0 = elem[i].vol;

				total_vol += elem[i].vol_0;

				// Keep a check on which element has the smallest volume
				if(elem[i].vol_0 < min_vol) {
					min_vol = elem[i].vol_0;
					min_vol_elem = i;
				}

				// Calc the mass of the element
				elem[i].mass = elem[i].vol_0 * elem[i].rho;
			}
			if(blob_state == FFEA_BLOB_IS_STATIC) {
				printf("\t\tBlob is static, so volume not defined within simulation.\n");
				printf("\t\tDefining Total rest volume of Blob to be 0 cubic Angstroms.\n");
				printf("\t\tAll elements have volume 0 cubic Angstroms.\n");	
				min_vol = 0.0;
 			} else {
				printf("\t\tTotal rest volume of Blob is %e cubic Angstroms.\n", total_vol * 1e30);
				printf("\t\tSmallest element (%i) has volume %e cubic Angstroms.\n", min_vol_elem, min_vol * 1e30);
			}

			// Calc the total mass of this Blob
			for(int i = 0; i < num_elements; i++) {
				mass += elem[i].mass;
			}
			printf("\t\tTotal mass of blob = %e\n", mass);
		}

		/*
 		*
 		*/
		int aggregate_forces_and_solve()
		{
			int n, m;

			// Aggregate the forces on each node by summing the contributions from each element.
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) private(n, m) schedule(guided)
			#endif
			for(n = 0; n < num_nodes; n++) {
				force[n].x = 0;
				force[n].y = 0;
				force[n].z = 0;
				for(m = 0; m < node[n].num_element_contributors; m++) {
					force[n].x += node[n].force_contributions[m]->x;
					force[n].y += node[n].force_contributions[m]->y;
					force[n].z += node[n].force_contributions[m]->z;
				}
			}

			// Aggregate surface forces onto nodes
			for(n = 0; n < num_surface_faces; n++) {
				for(int i = 0; i < 3; i++) {
					int sni = surface[n].n[i]->index;
					force[sni].x += surface[n].force[i].x;
					force[sni].y += surface[n].force[i].y;
					force[sni].z += surface[n].force[i].z;

		//			printf("force on %d from face %d = %e %e %e\n", sni, n, force[sni].x, force[sni].y, force[sni].z);
				}
			}
		//	printf("----\n\n");

			if(params->do_stokes == 1) {
				if(linear_solver != FFEA_NOMASS_CG_SOLVER) {
					#ifdef FFEA_PARALLEL_WITHIN_BLOB
						#pragma omp parallel for default(none) schedule(guided)
					#endif
					for(int i = 0; i < num_nodes; i++) {
						int thread_id = omp_get_thread_num();
						force[i].x -= node[i].vel.x * node[i].stokes_drag;
						force[i].y -= node[i].vel.y * node[i].stokes_drag;
						force[i].z -= node[i].vel.z * node[i].stokes_drag;
						if(params->calc_noise == 1) {
							force[i].x -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
							force[i].y -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
							force[i].z -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
						}
					}
				} else {
					if(params->calc_noise == 1) {
						#ifdef FFEA_PARALLEL_WITHIN_BLOB
							#pragma omp parallel for default(none) schedule(guided)
						#endif
						for(int i = 0; i < num_nodes; i++) {
							int thread_id = omp_get_thread_num();
							force[i].x -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
							force[i].y -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
							force[i].z -= RAND(-.5,.5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
						}
					}
				}
			}

		//	printf("2 num_nodes = %d\n", num_nodes);
		//	for(int i = 0; i < num_nodes; i++) {
		//		printf("%e %e %e\n", force[i].x, force[i].y, force[i].z);
		//	}


			// Set to zero any forces on the pinned nodes
			for(n = 0; n < num_pinned_nodes; n++) {
				int pn_index = pinned_nodes_list[n];
				force[pn_index].x = 0;
				force[pn_index].y = 0;
				force[pn_index].z = 0;
			}
			
			// Use the linear solver to solve for Mx = f where M is the Blob's mass matrix, 
			// or Kv = f where K is the viscosity matrix for the system
			// x/v is the (unknown) force solution and f is the force vector for the system.
			if(solver->solve(force) == FFEA_ERROR) {
				FFEA_ERROR_MESSG("Error reported by Solver.\n");
			}

			return FFEA_OK;
		}

		/*
 		 *
 		 */
		void euler_integrate()
		{
			int i;
			
			// Update the velocities and positions of all the nodes
			if(linear_solver == FFEA_NOMASS_CG_SOLVER) {
				#ifdef FFEA_PARALLEL_WITHIN_BLOB
					#pragma omp parallel for default(none) private(i) schedule(static)
				#endif
				for(i = 0; i < num_nodes; i++) {
					node[i].vel.x = force[i].x;
					node[i].vel.y = force[i].y;
					node[i].vel.z = force[i].z;

					node[i].pos.x += node[i].vel.x * params->dt;
					node[i].pos.y += node[i].vel.y * params->dt;
					node[i].pos.z += node[i].vel.z * params->dt;
				}

			} else {
				#ifdef FFEA_PARALLEL_WITHIN_BLOB
					#pragma omp parallel for default(none) private(i) schedule(static)
				#endif
				for(i = 0; i < num_nodes; i++) {
					node[i].vel.x += force[i].x * params->dt;
					node[i].vel.y += force[i].y * params->dt;
					node[i].vel.z += force[i].z * params->dt;

					node[i].pos.x += node[i].vel.x * params->dt;
					node[i].pos.y += node[i].vel.y * params->dt;
					node[i].pos.z += node[i].vel.z * params->dt;
				}
			}
		}

		/*
 		 *
		 */
		int calculate_node_element_connectivity()
		{
			int i, j;
			int *node_counter = NULL;

			// initialise num_element_contributors to zero
			for(i = 0; i < num_nodes; i++)
				node[i].num_element_contributors = 0;

			// count how many times each node is referenced in the list of elements
			for(i = 0; i < num_elements; i++)
				for(j = 0; j < NUM_NODES_QUADRATIC_TET; j++)
					node[elem[i].n[j]->index].num_element_contributors++;

			// allocate the contributions array for each node to the length given by num_element_contributors
			for(i = 0; i < num_nodes; i++) {
				node[i].force_contributions = new vector3 *[node[i].num_element_contributors];
				if(node[i].force_contributions == NULL) {
					FFEA_ERROR_MESSG("Failed to allocate memory for 'force_contributions' array (on node %d)\n", i);
				}
			}

			// create an array of counters keeping track of how full each contributions array is on each node
			node_counter = new int[num_nodes];
			if(node_counter == NULL) {
				FFEA_ERROR_MESSG("Failed to allocate node_counter array.\n");
			}

			// initialise the node_counter array to zero
			for(i = 0; i < num_nodes; i++)
				node_counter[i] = 0;

			// go back through the elements array and fill the contributions arrays with pointers to the
			// appropriate force contributions in the elements
			int node_index;
			for(i = 0; i < num_elements; i++)
				for(j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
					node_index = elem[i].n[j]->index;
					node[node_index].force_contributions[node_counter[node_index]] = &elem[i].node_force[j];
					node_counter[node_index]++;
				}

			// release node counter array
			delete[] node_counter;

			return FFEA_OK;
		}


		void build_mass_matrix()
		{
			// Calculate the Sparsity Pattern for the Mass matrix
			printf("\t\tCalculating sparsity pattern for 2nd order Mass matrix\n");
			SparsityPattern sparsity_pattern_mass_matrix;
			sparsity_pattern_mass_matrix.init(num_nodes);

			MassMatrixQuadratic M_alpha[num_elements];

			scalar *mem_loc;
			int ni_index, nj_index;
			for(int el = 0; el < num_elements; el++) {

				elem[el].construct_element_mass_matrix(&M_alpha[el]);

				for(int ni = 0; ni < 10; ni++) {
					for(int nj = 0; nj < 10; nj++) {

						ni_index = elem[el].n[ni]->index;
						nj_index = elem[el].n[nj]->index;

						mem_loc = M_alpha[el].get_M_alpha_mem_loc(ni, nj);

						sparsity_pattern_mass_matrix.register_contribution(ni_index, nj_index, mem_loc);
					}
				}
			}

			printf("\t\tBuilding sparsity pattern\n");

			// Use the sparsity patterns to create a fixed pattern sparse matrix
			M = sparsity_pattern_mass_matrix.create_sparse_matrix();

			// Build the mass matrix
			M->build();
		}
};

#endif
