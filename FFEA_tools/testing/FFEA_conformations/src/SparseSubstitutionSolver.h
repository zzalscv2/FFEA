#ifndef SPARSESUBSTITUTIONSOLVER_H_INCLUDED
#define SPARSESUBSTITUTIONSOLVER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"

#define INDEX(I, J) ((I) * num_rows + (J))

class SparseSubstitutionSolver: public Solver
{
	public:
		/* Constructor */
		SparseSubstitutionSolver()
		{
			// Initialise everything to zero
			num_rows = 0;
			L_key = NULL;
			U_key = NULL;
			inverse_diag = NULL;
			L = NULL;
			U = NULL;
		}

		/* Destructor */
		~SparseSubstitutionSolver()
		{
			num_rows = 0;
			delete[] L_key;
			delete[] U_key;
			delete[] inverse_diag;
			delete[] L;
			delete[] U;
			L_key = NULL;
			U_key = NULL;
			inverse_diag = NULL;
			L = NULL;
			U = NULL;
		}

		/* Builds the lower triangular Cholesky decomposed mass matrix */
		int init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list)
		{

			// Mass matrix will have as many rows as there are nodes in the mesh
			num_rows = num_nodes;

			printf("Allocating memory for mass matrix...\n");
			scalar *mass = new scalar[num_rows * num_rows];
			scalar *mass_LU = new scalar[num_rows * num_rows];
			printf("...done.\n");

			printf("Zeroing...\n");
			for(int i = 0; i < num_rows * num_rows; i++) {
				mass[i] = 0;
				mass_LU[i] = 0;
			}
			printf("...done\n");

			// Create a temporary lookup for checking if a node is 'pinned' or not.
			// if it is, then only a 1 on the diagonal corresponding to that node should
			// be placed (no off diagonal), effectively taking this node out of the equation
			// and therefore meaning the force on it should always be zero.
			int is_pinned[num_nodes];
			for(int i = 0; i < num_nodes; i++) {
				is_pinned[i] = 0;
			}
			for(int i = 0; i < num_pinned_nodes; i++) {
				is_pinned[pinned_nodes_list[i]] = 1;
			}

			// build the matrix
			printf("Building the mass matrix...\n");
			int ni, nj;
		        for(int n = 0; n < num_elements; n++) {
                		// add mass matrix for this element
		                for(int i = 0; i < 10; i++) {
		                        for(int j = 0; j < 10; j++) {
						if(i < 4 && j < 4) {
							ni = elem[n].n[i]->index;
							nj = elem[n].n[j]->index;
							if(is_pinned[ni] == 0 && is_pinned[nj] == 0) {
				                                if(i == j) {
			        	                                mass[ni * num_rows + nj] += .1 * elem[n].rho * elem[n].vol_0;
			                	                } else {
									mass[ni * num_rows + nj] += .05 * elem[n].rho * elem[n].vol_0;
								}
							} else {
								if(i == j) {
			        	                                mass[ni * num_rows + nj] = 1;
								}
							}
						} else {
							if(i == j) {
								ni = elem[n].n[i]->index;
								nj = elem[n].n[j]->index;
								mass[ni * num_rows + nj] = 1;
							}
						}
					}
				}
			}
			printf("...done\n");

			/* Perform cholesky decomposition on the calculated mass matrix, storing the result in mass_LU */
			printf("Performing Cholesky decomposition...\n");
			for(int k = 0; k < num_rows; k++) {
				for(int i = 0; i < k + 1; i++) {
					scalar sum_ij_kj = 0;
					for(int j = 0; j < i; j++) {
						sum_ij_kj += mass_LU[INDEX(i,j)] * mass_LU[INDEX(k,j)];
					}
					mass_LU[INDEX(k,i)] = (mass[INDEX(k,i)] - sum_ij_kj)/mass_LU[INDEX(i,i)];
				}

				scalar sum_kj_2 = 0;
				for(int j = 0; j < k; j++) {
					sum_kj_2 += mass_LU[INDEX(k, j)] * mass_LU[INDEX(k, j)];
				}
				mass_LU[INDEX(k, k)] = sqrt(mass[INDEX(k, k)] - sum_kj_2);
			}
			delete[] mass;
			printf("...done.\n");

			// Copy the lower matrix into upper
			for(int i = 1; i < num_rows; i++) {
				for(int j = 0; j < i; j++) {
					mass_LU[j * num_rows + i] = mass_LU[i * num_rows + j];
				}
			}

        		// Allocate and fill the 'inverse_diag' array
			inverse_diag = new scalar[num_rows];

			for(int i = 0; i < num_rows; i++)
				inverse_diag[i] = 1.0/mass_LU[i * num_rows + i];

			// Allocate the 'key' arrays
			L_key = new int[num_rows];
			U_key = new int[num_rows];

			// Build the lower triangular matrix key
			int total_L = 0;
			for(int i = 0; i < num_rows; i++) {
				for(int j = num_rows-1; j >= i; j--) {
					if(mass_LU[i * num_rows + j] != 0) {
						total_L += (j - i);
						L_key[i] = (j - i);
						break;
					}
				}
			}

			// Build the upper triangular matrix key
			total_entries_in_U = 0;
			for(int i = 0; i < num_rows; i++)
				for(int j = 0; j <= i; j++)
					if(mass_LU[i * num_rows + j] != 0) {
						total_entries_in_U += (i - j);
						U_key[i] = (i - j);
						break;
					}

			// Allocate the off-diagonal entry arrays
			L = new scalar[total_L];
			U = new scalar[total_entries_in_U];

			// Fill up the L and U sparse triangular matrices
			int off_diag_data_index = 0;
			for(int i = 0; i < num_rows; i++)
				for(int j = i + 1; j <= i + L_key[i]; j++) {
					L[off_diag_data_index] = mass_LU[i * num_rows + j];
					off_diag_data_index++;
				}
			off_diag_data_index = 0;
			for(int i = 0; i < num_rows; i++)
				for(int j = i - U_key[i]; j < i; j++) {
					U[off_diag_data_index] = mass_LU[i * num_rows + j];
					off_diag_data_index++;
				}

			delete[] mass_LU;

			return FFEA_OK;
		}

		/*
		 * Solves the equation Ax = b for the unknown vector x, for
		 * 3 right-hand-sides (b vectors) at once, using forward and backward substitution. This function
		 * considers the case where A is sparse symmetric positive-definite and therefore the solver must
		 * have constructed a variable band matrix containing the result of a Cholesky decomposition on
		 * the original matrix A.
		 *
		 * This implementation is memory efficient as it modifies the vector x in place (no temporary vectors are needed).
		 */
		int solve(vector3 *x)
		{
			int i, j, index;

			// Forward substitution step Ly = b :
			index = 0;
			for(i = 0; i < num_rows; i++) {

				x[i].x *= inverse_diag[i];
				x[i].y *= inverse_diag[i];
				x[i].z *= inverse_diag[i];

				for(j = 0; j < L_key[i]; j++) {
					x[i + j + 1].x -= x[i].x * L[j + index];
					x[i + j + 1].y -= x[i].y * L[j + index];
					x[i + j + 1].z -= x[i].z * L[j + index];
				}

				index += L_key[i];
			}

			// Backward substitution step Ux = y :
			index = total_entries_in_U - 1;
			for(i = num_rows - 1; i >= 0; i--) {

				x[i].x *= inverse_diag[i];
				x[i].y *= inverse_diag[i];
				x[i].z *= inverse_diag[i];

				for(j = 0; j < U_key[i]; j++) {
					x[i - j - 1].x -= x[i].x * U[index - j];
					x[i - j - 1].y -= x[i].y * U[index - j];
					x[i - j - 1].z -= x[i].z * U[index - j];
				}

				index -= U_key[i];
			}

			return FFEA_OK;
		}

		void apply_matrix(scalar *in, scalar *result)
		{
		}

	private:

		/* Number of rows in this cholesky variable band matrix */
		int num_rows;

		/*
		 * Array of size num_rows. Each element contains the index of the last non-zero
		 * entry in the corresponding column. This is the end of the "band" of entries stored
		 * in that column, that starts from the diagonal (travelling downwards or upwards for
		 * the L and U matrices respectively).
		 */
		int *L_key, *U_key;

		/* Stores the number of entries stored in the upper triangle */
		int total_entries_in_U;

		/* Stores the inverse of the diagonal elements (need for LU solving) */
		scalar *inverse_diag;

		/*
		 * Stores all the data from the bands of the band matrix (note that this may contain
		 * some zeroes, but if it is a proper sparse band matrix there should be very few of these)
		 * NOT INCLUDING THE DIAGONAL: only the inverse of the diagonal.
		 *
		 * The off-diagonal data for the lower and upper triangles are stored in L and U respectively.
		 */
		scalar *L, *U;
};

#endif
