/*
 * tetra_element_linear.hpp
 */

#ifndef TETRA_ELEMENT_LINEAR_H_INCLUDED
#define TETRA_ELEMENT_LINEAR_H_INCLUDED

#include <stddef.h>
#include "mat_vec_fns.h"
#include "SimulationParams.h"
#include "FFEA_return_codes.h"
#include "MersenneTwister.h"
#include "SecondOrderFunctions.h"
#include "PoissonMatrixQuadratic.h"
#include "MassMatrixQuadratic.h"

class Blob;

#define NUM_NODES_QUADRATIC_TET 10

/*
 * A preprocessor macro generating code that calculates a submatrix of the bulk viscous
 * matrix. 'I' specifies what position in the matrix 'V' this 4x4 submatrix should be
 * inserted to (position I,I in the matrix).
 *
 * N.B. This macro only deals with submatrices *on the diagonal* of the bulk viscous matrix.
 */
#define BULK_VISCOUS_SUBMATRIX_DIAG(V,I) \
	V[I + 0][I + 0] = dpsi[I + 0] * dpsi[I + 0]	* (A + B); \
	V[I + 1][I + 1] = dpsi[I + 1] * dpsi[I + 1]	* (A + B); \
	V[I + 2][I + 2] = dpsi[I + 2] * dpsi[I + 2]	* (A + B); \
	V[I + 3][I + 3] = dpsi[I + 3] * dpsi[I + 3]	* (A + B); \
	V[I + 0][I + 1] = dpsi[I + 1] * dpsi[I + 0]	* (A + B); \
	V[I + 0][I + 2] = dpsi[I + 2] * dpsi[I + 0]	* (A + B); \
	V[I + 0][I + 3] = dpsi[I + 3] * dpsi[I + 0]	* (A + B); \
	V[I + 1][I + 2] = dpsi[I + 1] * dpsi[I + 2]	* (A + B); \
	V[I + 1][I + 3] = dpsi[I + 1] * dpsi[I + 3]	* (A + B); \
	V[I + 2][I + 3] = dpsi[I + 2] * dpsi[I + 3]	* (A + B);

/*
 * A preprocessor macro generating code that calculates a submatrix of the bulk viscous
 * matrix. 'I' and 'J' specify what position in the matrix 'V' this 4x4 submatrix should
 * be inserted to (position I,J in the matrix).
 *
 * N.B. This code will calculate submatrices anywhere in the bulk viscous matrix,
 * but the BULK_VISCOUS_SUBMATRIX_DIAG macro above provides simpler more efficient
 * code for the case in which the submatrix lies on the diagonal of the bulk viscous
 * matrix.
 */
#define BULK_VISCOUS_SUBMATRIX_OFFDIAG(V, I, J) \
	K[0][0] = dpsi[I] * dpsi[J]; \
	K[0][1] = dpsi[J + 1] * dpsi[I]; \
	K[0][2] = dpsi[J + 2] * dpsi[I]; \
	K[0][3] = dpsi[J + 3] * dpsi[I]; \
	\
	K[1][0] = dpsi[I + 1] * dpsi[J]; \
	K[1][1] = dpsi[I + 1] * dpsi[J + 1]; \
	K[1][2] = dpsi[I + 1] * dpsi[J + 2]; \
	K[1][3] = dpsi[I + 1] * dpsi[J + 3]; \
	\
	K[2][0] = dpsi[I + 2] * dpsi[J]; \
	K[2][1] = dpsi[I + 2] * dpsi[J + 1]; \
	K[2][2] = dpsi[I + 2] * dpsi[J + 2]; \
	K[2][3] = dpsi[I + 2] * dpsi[J + 3]; \
	\
	K[3][0] = dpsi[I + 3] * dpsi[J]; \
	K[3][1] = dpsi[I + 3] * dpsi[J + 1]; \
	K[3][2] = dpsi[I + 3] * dpsi[J + 2]; \
	K[3][3] = dpsi[I + 3] * dpsi[J + 3]; \
	\
	V[I + 0][J + 0] = (B + A) * K[0][0]; \
	V[I + 0][J + 1] = B * K[0][1] + A * K[1][0]; \
	V[I + 0][J + 2] = B * K[0][2] + A * K[2][0]; \
	V[I + 0][J + 3] = B * K[0][3] + A * K[3][0]; \
	\
	V[I + 1][J + 0] = B * K[1][0] + A * K[0][1]; \
	V[I + 1][J + 1] = (B + A) * K[1][1]; \
	V[I + 1][J + 2] = B * K[1][2] + A * K[2][1]; \
	V[I + 1][J + 3] = B * K[1][3] + A * K[3][1]; \
	\
	V[I + 2][J + 0] = B * K[2][0] + A * K[0][2]; \
	V[I + 2][J + 1] = B * K[2][1] + A * K[1][2]; \
	V[I + 2][J + 2] = (B + A) * K[2][2]; \
	V[I + 2][J + 3] = B * K[2][3] + A * K[3][2]; \
	\
	V[I + 3][J + 0] = B * K[3][0] + A * K[0][3]; \
	V[I + 3][J + 1] = B * K[3][1] + A * K[1][3]; \
	V[I + 3][J + 2] = B * K[3][2] + A * K[2][3]; \
	V[I + 3][J + 3] = (B + A) * K[3][3]; \

#define RAND(A, B) ((A) + ((B)-(A))*(rng[thread_id].rand()))

#define DPSI1_DX 0
#define DPSI2_DX 1
#define DPSI3_DX 2
#define DPSI4_DX 3
#define DPSI1_DY 4
#define DPSI2_DY 5
#define DPSI3_DY 6
#define DPSI4_DY 7
#define DPSI1_DZ 8
#define DPSI2_DZ 9
#define DPSI3_DZ 10
#define DPSI4_DZ 11

/*
 * A 10-point "quadratic" tetrahedron
 */
class tetra_element_linear
{

	public:

		tetra_element_linear()
		{
			rho = 0;
			A = 0;
			B = 0;
			G = 0;
			E = 0;
			dielectric = 0;
			for(int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				n[i] = NULL;
			}

//			node_phi[0] = 0; node_phi[1] = 0; node_phi[2] = 0; node_phi[3] = 0;
			vol_0 = 0;
			vol = 0;
			mat3_set_zero(F_ij);
			internal_stress_mag = 0;
			mat3_set_zero(J_inv_0);
			mat12_set_zero(viscosity_matrix);
			zero_force();
			last_det = 0;
		}

		/*
		 * Get the memory location of the specified element of K_alpha
		 */
		scalar * get_K_alpha_element_mem_loc(int ni, int nj)
		{
			return K_alpha.get_K_alpha_mem_loc(ni, nj);
		}

		/* Calc the diffusion matrix for this element */
		void calculate_K_alpha()
		{
			// Build the poisson diffusion matrix corresponding to this element
			K_alpha.build(n, dielectric);
		}

		void construct_element_mass_matrix(MassMatrixQuadratic *M_alpha)
		{
			// Build the element mass matrix corresponding to this element
			M_alpha->build(n);
		}

		void add_K_alpha(scalar *K, int num_nodes)
		{
			for(int i = 0; i < 10; i++) {
				for(int j = 0; j < 10; j++) {
					K[n[i]->index * num_nodes + n[j]->index] += K_alpha.get_K_alpha_value(i, j);
				}
			}
		}

		/* Returns the gradient of the potential at the given (s,t,u) position in the element */
		void get_grad_phi_at_stu(vector3 *grad_phi, scalar s, scalar t, scalar u)
		{
			vector3 grad_psi[NUM_NODES_QUADRATIC_TET];

			SecondOrderFunctions::abcd J_coeff[3][3];
			SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);
			scalar J_inv[9];
			scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, s, t, u, J_inv);
			SecondOrderFunctions::calc_grad_psi(grad_psi, s, t, u, J_inv);

			grad_phi->x = 0;
			grad_phi->y = 0;
			grad_phi->z = 0;
			for(int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				grad_phi->x += grad_psi[i].x * n[i]->phi;
				grad_phi->y += grad_psi[i].y * n[i]->phi;
				grad_phi->z += grad_psi[i].z * n[i]->phi;
			}

			grad_phi->x *= det_J;
			grad_phi->y *= det_J;
			grad_phi->z *= det_J;
		}

		/* Calculates the force on each node of the element due to the electrostatic potential gradient there */
		void calculate_electrostatic_forces()
		{
			struct tetrahedron_gauss_point gauss_points[NUM_TET_GAUSS_QUAD_POINTS] =
			{
				// Weight, eta1, eta2, eta3, eta4
				{0.317460317460317450e-2, {0.5, 0.5, 0.0, 0.0}},
				{0.317460317460317450e-2, {0.5, 0.0, 0.5, 0.0}},
				{0.317460317460317450e-2, {0.5, 0.0, 0.0, 0.5}},
				{0.317460317460317450e-2, {0.0, 0.5, 0.5, 0.0}},
				{0.317460317460317450e-2, {0.0, 0.0, 0.5, 0.5}},
				{0.317460317460317450e-2, {0.0, 0.5, 0.0, 0.5}},
				{0.147649707904967828e-1, {0.100526765225204467, 0.100526765225204467, 0.100526765225204467, 0.698419704324386603}},
				{0.147649707904967828e-1, {0.100526765225204467, 0.100526765225204467, 0.698419704324386603, 0.100526765225204467}},
				{0.147649707904967828e-1, {0.100526765225204467, 0.698419704324386603, 0.100526765225204467, 0.100526765225204467}},
				{0.147649707904967828e-1, {0.698419704324386603, 0.100526765225204467, 0.100526765225204467, 0.100526765225204467}},
				{0.221397911142651221e-1, {0.314372873493192195, 0.314372873493192195, 0.314372873493192195, 0.568813795204234229e-1}},
				{0.221397911142651221e-1, {0.314372873493192195, 0.314372873493192195, 0.568813795204234229e-1, 0.314372873493192195}},
				{0.221397911142651221e-1, {0.314372873493192195, 0.568813795204234229e-1, 0.314372873493192195, 0.314372873493192195}},
				{0.221397911142651221e-1, {0.568813795204234229e-1, 0.314372873493192195, 0.314372873493192195, 0.314372873493192195}}
			};

			SecondOrderFunctions::abcd J_coeff[3][3];
			vector3 grad_psi[NUM_NODES_QUADRATIC_TET];
			scalar psi[NUM_NODES_QUADRATIC_TET];
			vector3 force[NUM_NODES_QUADRATIC_TET];
			vector3 grad_phi_here;

			for(int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				force[i].x = 0;
				force[i].y = 0;
				force[i].z = 0;
			}

			SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);

			scalar J_inv[9];

			for(int i = 0; i < NUM_TET_GAUSS_QUAD_POINTS; i++) {
				scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
				SecondOrderFunctions::calc_grad_psi(grad_psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
				SecondOrderFunctions::calc_psi(psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2]);

				grad_phi_here.x = 0;
				grad_phi_here.y = 0;
				grad_phi_here.z = 0;
				for(int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
					grad_phi_here.x += grad_psi[j].x * n[j]->phi * det_J;
					grad_phi_here.y += grad_psi[j].y * n[j]->phi * det_J;
					grad_phi_here.z += grad_psi[j].z * n[j]->phi * det_J;
//					printf("%d >> %e (%e %e %e)\n", i, det_J, grad_psi[j].x, grad_psi[j].y, grad_psi[j].z);
				}

//				printf("%d :: (%e, %e, %e)\n", i, grad_phi_here.x, grad_phi_here.y, grad_phi_here.z);

				for(int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
//					force[j].x -= gauss_points[i].W * grad_phi_here.x * psi[j] * n[j]->rho;
//					force[j].y -= gauss_points[i].W * grad_phi_here.y * psi[j] * n[j]->rho;
//					force[j].z -= gauss_points[i].W * grad_phi_here.z * psi[j] * n[j]->rho;
//					printf("%d == %e\n", i, n[j]->rho);
					force[j].x -= gauss_points[i].W * grad_phi_here.x * n[j]->rho;
					force[j].y -= gauss_points[i].W * grad_phi_here.y * n[j]->rho;
					force[j].z -= gauss_points[i].W * grad_phi_here.z * n[j]->rho;
				}
			}

			for(int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				add_force_to_node(i, &force[i]);
			}

		}

		/* Calculates the Jacobian matrix for this element */
		void calculate_jacobian(matrix3 J)
		{
			J[0][0] = n[1]->pos.x - n[0]->pos.x;
			J[0][1] = n[1]->pos.y - n[0]->pos.y;
			J[0][2] = n[1]->pos.z - n[0]->pos.z;

			J[1][0] = n[2]->pos.x - n[0]->pos.x;
			J[1][1] = n[2]->pos.y - n[0]->pos.y;
			J[1][2] = n[2]->pos.z - n[0]->pos.z;

			J[2][0] = n[3]->pos.x - n[0]->pos.x;
			J[2][1] = n[3]->pos.y - n[0]->pos.y;
			J[2][2] = n[3]->pos.z - n[0]->pos.z;
		}


		/*
		 * Inverts the given jacobian matrix J, using this to calculate the derivatives
		 * of the shape functions which are stored in dpsi. This is an array
		 * of all 12 derivatives, in the following order:
		 *		dpsi    =	[ d(psi)_1/dx ]
		 *				[ d(psi)_2/dx ]
		 *				[ d(psi)_3/dx ]
		 *				[ d(psi)_4/dx ]
		 *				[ d(psi)_1/dy ]
		 *				[ d(psi)_2/dy ]
		 *				[     ...     ]
		 *				[ d(psi)_4/dz ]
		 *
		 * Function also (as a by-product of the inversion) calculates the volume of the
		 * element whose jacobian this is, which is stored in 'vol'.
		 */
		int calc_shape_function_derivatives_and_volume(matrix3 J)
		{
			scalar det;

			// Calculate shape function derivs from inverse jacobian directly into dpsi[]

/*
			dpsi[1] = J[2][2]*J[1][1] - J[2][1]*J[1][2];
			dpsi[2] = J[2][1]*J[0][2] - J[2][2]*J[0][1];
			dpsi[3] = J[1][2]*J[0][1] - J[1][1]*J[0][2];

			dpsi[5] = J[2][0]*J[1][2] - J[2][2]*J[1][0];
			dpsi[6] = J[2][2]*J[0][0] - J[2][0]*J[0][2];
			dpsi[7] = J[1][0]*J[0][2] - J[1][2]*J[0][0];

			dpsi[9] = J[2][1]*J[1][0] - J[2][0]*J[1][1];
			dpsi[10] = J[2][0]*J[0][1] - J[2][1]*J[0][0];
			dpsi[11] = J[1][1]*J[0][0] - J[1][0]*J[0][1];
*/
			dpsi[DPSI2_DX] = J[2][2]*J[1][1] - J[2][1]*J[1][2];
			dpsi[DPSI3_DX] = J[2][1]*J[0][2] - J[2][2]*J[0][1];
			dpsi[DPSI4_DX] = J[1][2]*J[0][1] - J[1][1]*J[0][2];

			dpsi[DPSI2_DY] = J[2][0]*J[1][2] - J[2][2]*J[1][0];
			dpsi[DPSI3_DY] = J[2][2]*J[0][0] - J[2][0]*J[0][2];
			dpsi[DPSI4_DY] = J[1][0]*J[0][2] - J[1][2]*J[0][0];

			dpsi[DPSI2_DZ] = J[2][1]*J[1][0] - J[2][0]*J[1][1];
			dpsi[DPSI3_DZ] = J[2][0]*J[0][1] - J[2][1]*J[0][0];
			dpsi[DPSI4_DZ] = J[1][1]*J[0][0] - J[1][0]*J[0][1];

			// Determinant
			det = J[0][0] * dpsi[DPSI2_DX] + J[1][0] * dpsi[DPSI3_DX] + J[2][0] * dpsi[DPSI4_DX];

			// Check if element has inverted itself (determinant changed sign)
			if(last_det * det < 0) {
				return FFEA_ERROR;
			}
			last_det = det;

			// Calculate volume of element from the determinant
			// and the volume of the "parent" element in local coords
			// e.g 1/6 for a right angled tetrahedron
			vol = (1.0/6.0) * fabs(det);

			// Divide by determinant to complete matrix inversion
			det = 1.0/det;
			dpsi[1] *= det;	dpsi[2] *= det;	dpsi[3] *= det;
			dpsi[5] *= det;	dpsi[6] *= det;	dpsi[7] *= det;
			dpsi[9] *= det;	dpsi[10] *= det; dpsi[11] *= det;

			// calculate the derivatives of shape function psi_1 w.r.t. x, y and z
/*
			dpsi[0] = -(dpsi[1] + dpsi[2] + dpsi[3]);
			dpsi[4] = -(dpsi[5] + dpsi[6] + dpsi[7]);
			dpsi[8] = -(dpsi[9] + dpsi[10] + dpsi[11]);
*/
			dpsi[DPSI1_DX] = -(dpsi[DPSI2_DX] + dpsi[DPSI3_DX] + dpsi[DPSI4_DX]);
			dpsi[DPSI1_DY] = -(dpsi[DPSI2_DY] + dpsi[DPSI3_DY] + dpsi[DPSI4_DY]);
			dpsi[DPSI1_DZ] = -(dpsi[DPSI2_DZ] + dpsi[DPSI3_DZ] + dpsi[DPSI4_DZ]);

			return FFEA_OK;
		}

		/*
		 * Builds the viscosity matrix from the shape function derivatives, the shear and bulk
		 * viscosity constants, and the element volume
		 */
		void create_viscosity_matrix()
		{
			int i, j;
			matrix4 K;

			// Construct submatrices on the diagonal
			BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix,0)
			BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix,4)
			BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix,8)

			// Construct submatrices off diagonal
			BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix,0,4)
			BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix,0,8)
			BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix,4,8)

			// Create the diffusion matrix and add it to the viscosity matrix in 3 upper
			// triangular blocks along the diagonal
			calc_del2_matrix();
			add_diffusion_matrix(viscosity_matrix);

			// Multiply all these values (currently only in upper half of matrix) by the element volume
			for(i = 0; i < 12; i++)
				for(j = 0; j <= i; j++)
					viscosity_matrix[j][i] *= vol;

			// Viscosity matrix is symmetric, so no need to recalculate entries.
			// Simply 'mirror image' the matrix back into itself
			for(i = 1; i < 12; i++)
				for(j = 0; j < i; j++)
					viscosity_matrix[i][j] = viscosity_matrix[j][i];
		}

		/*
		 *
		 */
		void add_shear_elastic_stress(matrix3 J, matrix3 stress)
		{
			// Reset gradient deformation to zero
			mat3_set_zero(F_ij);

			// F_ij transpose is given by the current jacobian times the rest state jacobian inverse.
			mat3_mult_both_transposed(J, J_inv_0, F_ij);

			// Now multiply F_ij with its own transpose to get F_ik_kjT
			mat3_mult_transpose(F_ij, F_ij, stress);

			// Scale matrix by the appropriate factors
			mat3_scale(stress, (G * vol_0/vol));

			// Subtract Identity to make deformation tensor (almost) traceless
			stress[0][0] -= G;
			stress[1][1] -= G;
			stress[2][2] -= G;
		}

		/*
		 *
		 */
		void add_bulk_elastic_stress(matrix3 stress)
		{
			//scalar c = (E + G/3.0) * ((vol - vol_0)/vol_0); //Old
			scalar c_2 = E - G * 2.0/3.0;
			scalar c = G * (1.0 - (vol_0/vol)) + 0.5 * c_2 * ((vol/vol_0) - (vol_0/vol));
			stress[0][0] += c;
			stress[1][1] += c;
			stress[2][2] += c;
		}

		/*
		 * Given the shape function derivatives, the element volume and a random number generator, this
		 * function calculates the fluctuating stress tensor, generating a stochastic change in the
		 * nodal velocities for the element under consideration. This function will add its contribution
		 * to the given 12-vector du.
		 *
		 */
		void add_fluctuating_stress(SimulationParams *params, MTRand rng[], matrix3 stress, int thread_id)
		{
			scalar c = sqrt((24 * params->kT)/(vol * params->dt));

			// Bulk fluctuation term
			scalar bf = sqrt(B) * RAND(-.5,.5);

			// Diagonal terms
			stress[0][0] += c * (sqrt(2 * A) * RAND(-.5,.5) + bf);
			stress[1][1] += c * (sqrt(2 * A) * RAND(-.5,.5) + bf);
			stress[2][2] += c * (sqrt(2 * A) * RAND(-.5,.5) + bf);

			// Off diagonal terms (note that stress matrix is symmetric)
			stress[0][1] += c * sqrt(A) * RAND(-.5,.5);
			stress[0][2] += c * sqrt(A) * RAND(-.5,.5);
			stress[1][2] += c * sqrt(A) * RAND(-.5,.5);

			stress[1][0] = stress[0][1];
			stress[2][0] = stress[0][2];
			stress[2][1] = stress[1][2];
		}

		/*
		 * Applies the given stress tensor to the shape function derivatives to get the contribution to du
		 */
		void apply_stress_tensor(matrix3 stress, vector12 du)
		{
			int i;
			for(i = 0; i < 3; i++) {
				du[4 * i] += vol * (dpsi[0] * stress[i][0] + dpsi[4] * stress[i][1] + dpsi[8] * stress[i][2]);
				du[4 * i + 1] += vol * (dpsi[1] * stress[i][0] + dpsi[5] * stress[i][1] + dpsi[9] * stress[i][2]);
				du[4 * i + 2] += vol * (dpsi[2] * stress[i][0] + dpsi[6] * stress[i][1] + dpsi[10] * stress[i][2]);
				du[4 * i + 3] += vol * (dpsi[3] * stress[i][0] + dpsi[7] * stress[i][1] + dpsi[11] * stress[i][2]);
			}
		}

		/*
		 * Sets the given 12-vector to the velocities of this element's four nodes,
		 */
		void get_element_velocity_vector(vector12 v)
		{
			v[0] = n[0]->vel.x;
			v[1] = n[1]->vel.x;
			v[2] = n[2]->vel.x;
			v[3] = n[3]->vel.x;

			v[4] = n[0]->vel.y;
			v[5] = n[1]->vel.y;
			v[6] = n[2]->vel.y;
			v[7] = n[3]->vel.y;

			v[8] = n[0]->vel.z;
			v[9] = n[1]->vel.z;
			v[10] = n[2]->vel.z;
			v[11] = n[3]->vel.z;
		}

		/*
		 * Add this element's nodal forces to those given in the force 12-vector
		 */
		void add_element_force_vector(vector12 force)
		{
			node_force[0].x -= force[0];
			node_force[1].x -= force[1];
			node_force[2].x -= force[2];
			node_force[3].x -= force[3];

			node_force[0].y -= force[4];
			node_force[1].y -= force[5];
			node_force[2].y -= force[6];
			node_force[3].y -= force[7];

			node_force[0].z -= force[8];
			node_force[1].z -= force[9];
			node_force[2].z -= force[10];
			node_force[3].z -= force[11];
		}

		/* Add given force to the specified node of this element */
		void add_force_to_node(int i, vector3 *f)
		{
			node_force[i].x += f->x;
			node_force[i].y += f->y;
			node_force[i].z += f->z;
		}

		/* A roundabout and inefficient way of working out what node (from 0 to 9) this index corresponds to */
		int what_node_is_this(int index)
		{
			int i;
			for(i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				if(n[i]->index == index) {
					return i;
				}
			}

			FFEA_ERROR_MESSG("Specified node index does not belong to this element\n")
		}

		void print()
		{
			int i;
			for(i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				printf("Node %d:\n", i);
				n[i]->print();
				printf("node_force: %e %e %e\n", node_force[i].x, node_force[i].y, node_force[i].z);
				printf("volume: %e\n", vol);
			}

		}

		/*
		 * Applies the mass matrix (for a linear tetrahedral element of density rho and equilibrium volume vol_0)
		 * to the force vector du to get the correct distribution of force between nodes.
		 */
		void apply_element_mass_matrix(vector12 du)
		{
			int i, j, k;

			// mass matrix
			const matrix4 M =      {{.1, .05, .05, .05},
						{.05, .1, .05, .05},
						{.05, .05, .1, .05},
						{.05, .05, .05, .1}};

			vector12 temp_du = {0,0,0,0,0,0,0,0,0,0,0,0};

			// Apply mass matrix to each of the three subvectors
			for(k = 0; k < 12; k += 4)
				for(i = 0; i < 4; i++)
					for(j = 0; j < 4; j++)
						temp_du[i + k] += M[i][j] * du[j + k];

			for(i = 0; i < 12; i++) du[i] = temp_du[i] * (vol_0 * rho);
		}

		void volume_coord_to_xyz(scalar eta0, scalar eta1, scalar eta2, scalar eta3, vector3 *r)
		{
			r->x = eta0 * n[0]->pos.x + eta1 * n[1]->pos.x + eta2 * n[2]->pos.x + eta3 * n[3]->pos.x;
			r->y = eta0 * n[0]->pos.y + eta1 * n[1]->pos.y + eta2 * n[2]->pos.y + eta3 * n[3]->pos.y;
			r->z = eta0 * n[0]->pos.z + eta1 * n[1]->pos.z + eta2 * n[2]->pos.z + eta3 * n[3]->pos.z;
		}

		void zero_force()
		{
			for(int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
				vector3_set_zero(&node_force[i]);
			}
		}

		void linearise_element()
		{

			n[4]->pos.x = .5 * (n[0]->pos.x + n[1]->pos.x);
			n[4]->pos.y = .5 * (n[0]->pos.y + n[1]->pos.y);
			n[4]->pos.z = .5 * (n[0]->pos.z + n[1]->pos.z);

			n[5]->pos.x = .5 * (n[0]->pos.x + n[2]->pos.x);
			n[5]->pos.y = .5 * (n[0]->pos.y + n[2]->pos.y);
			n[5]->pos.z = .5 * (n[0]->pos.z + n[2]->pos.z);

			n[6]->pos.x = .5 * (n[0]->pos.x + n[3]->pos.x);
			n[6]->pos.y = .5 * (n[0]->pos.y + n[3]->pos.y);
			n[6]->pos.z = .5 * (n[0]->pos.z + n[3]->pos.z);

			n[7]->pos.x = .5 * (n[1]->pos.x + n[2]->pos.x);
			n[7]->pos.y = .5 * (n[1]->pos.y + n[2]->pos.y);
			n[7]->pos.z = .5 * (n[1]->pos.z + n[2]->pos.z);

			n[8]->pos.x = .5 * (n[1]->pos.x + n[3]->pos.x);
			n[8]->pos.y = .5 * (n[1]->pos.y + n[3]->pos.y);
			n[8]->pos.z = .5 * (n[1]->pos.z + n[3]->pos.z);

			n[9]->pos.x = .5 * (n[2]->pos.x + n[3]->pos.x);
			n[9]->pos.y = .5 * (n[2]->pos.y + n[3]->pos.y);
			n[9]->pos.z = .5 * (n[2]->pos.z + n[3]->pos.z);
		}

		void calc_centroid()
		{
			centroid.x = .25 * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x + n[3]->pos.x);
			centroid.y = .25 * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y + n[3]->pos.y);
			centroid.z = .25 * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z + n[3]->pos.z);
		}

		/* Properties of this element */
		scalar rho;     /* density */
		scalar A;       /* shear viscosity */
		scalar B;       /* second coefficient of viscosity */
		scalar G;       /* shear modulus */
		scalar E;       /* bulk modulus */
		scalar dielectric;
		scalar mass;

		/* A quadratic tetrahedron has 10 nodes. Keep pointers to the actual memory location of the nodes. */
		mesh_node *n[NUM_NODES_QUADRATIC_TET];

		/* The 12-vector containing the shape function derivatives for this element */
		vector12 dpsi;

		/*
	 	 * The del2 matrix for this element (in upper triangular form, since it's symmetric).
		 * Used in constructing the diffusion matrix and the poisson matrix.
		 */
		upper_triangular_matrix4 del2;

		PoissonMatrixQuadratic K_alpha;

		/* Store the contribution from this element to the force on each of its four nodes */
		vector3 node_force[NUM_NODES_QUADRATIC_TET];

		/* The rest volume of this element */
		scalar vol_0;

		/* The current volume of this element */
		scalar vol;

		/* The gradient deformation tensor for this element (needed for potential energy calculation) */
		matrix3 F_ij;

		/* The double contraction of the internal stress tensor, including elastic and thermal stresses */
		scalar internal_stress_mag;

		/* The inverse jacobian of this element at rest */
		matrix3 J_inv_0;

		/* Viscosity Matrix for the internal forces of this element */
		matrix12 viscosity_matrix;

		/* Blob this element is a part of */
		Blob *daddy_blob;

		/* Index of this element in the parent Blob */
		int index;

		vector3 centroid;

	private:

		/*
		 * The last determinant of this element's transformation (used to work out whether it has inverted itself)
		 */
		scalar last_det;

		/*
		 * Creates del2 matrix from the shape function derivatives
		 */
		void calc_del2_matrix()
		{
			del2.u00 = (dpsi[0]*dpsi[0] + dpsi[4]*dpsi[4] + dpsi[8]*dpsi[8]);
			del2.u01 = (dpsi[0]*dpsi[1] + dpsi[4]*dpsi[5] + dpsi[8]*dpsi[9]);
			del2.u02 = (dpsi[0]*dpsi[2] + dpsi[4]*dpsi[6] + dpsi[8]*dpsi[10]);
			del2.u03 = (dpsi[0]*dpsi[3] + dpsi[4]*dpsi[7] + dpsi[8]*dpsi[11]);

			del2.u11 = (dpsi[1]*dpsi[1] + dpsi[5]*dpsi[5] + dpsi[9]*dpsi[9]);
			del2.u12 = (dpsi[1]*dpsi[2] + dpsi[5]*dpsi[6] + dpsi[9]*dpsi[10]);
			del2.u13 = (dpsi[1]*dpsi[3] + dpsi[5]*dpsi[7] + dpsi[9]*dpsi[11]);

			del2.u22 = (dpsi[2]*dpsi[2] + dpsi[6]*dpsi[6] + dpsi[10]*dpsi[10]);
			del2.u23 = (dpsi[2]*dpsi[3] + dpsi[6]*dpsi[7] + dpsi[10]*dpsi[11]);

			del2.u33 = (dpsi[3]*dpsi[3] + dpsi[7]*dpsi[7] + dpsi[11]*dpsi[11]);

//			printf("del2:\n");
//			printf("%e %e %e %e\n", del2.u00, del2.u01, del2.u02, del2.u03);
//			printf("%e %e %e %e\n", del2.u01, del2.u11, del2.u12, del2.u13);
//			printf("%e %e %e %e\n", del2.u02, del2.u12, del2.u22, del2.u23);
//			printf("%e %e %e %e\n", del2.u03, del2.u13, del2.u23, del2.u33);
		}

		void add_diffusion_matrix(matrix12 V)
		{
			int i;

			// Drop each upper triangle of this diffusion matrix along the block diagonal
			// of the viscosity matrix
			for(i = 0; i < 3; i++) {
				V[i*4+0][i*4+0]+=del2.u00 * A;  V[i*4+0][i*4+1]+=del2.u01 * A;  V[i*4+0][i*4+2]+=del2.u02 * A;  V[i*4+0][i*4+3]+=del2.u03 * A;
								V[i*4+1][i*4+1]+=del2.u11 * A;  V[i*4+1][i*4+2]+=del2.u12 * A;  V[i*4+1][i*4+3]+=del2.u13 * A;
												V[i*4+2][i*4+2]+=del2.u22 * A;  V[i*4+2][i*4+3]+=del2.u23 * A;
																V[i*4+3][i*4+3]+=del2.u33 * A;
			}
		}

		struct tetrahedron_gauss_point
		{
			scalar W;
			scalar eta[4];
		};
};

#endif
