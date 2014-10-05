#ifndef MASSMATRIXQUADRATIC_H_INCLUDED
#define MASSMATRIXQUADRATIC_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_10X10 55

#define NUM_SHAPE_FUNCTIONS 10

#define NUM_TET_GAUSS_QUAD_POINTS 14

#include "SecondOrderFunctions.h"


class MassMatrixQuadratic
{
	public:
		MassMatrixQuadratic()
		{
			zero();
		}

		scalar * get_M_alpha_mem_loc(int i, int j)
		{
			// Poisson matrix is symmetric, so convert any request for an upper triangular element into its
			// corresponding (equivalent) lower triangular element
			if(i < j) {
				int temp = i;
				i = j;
				j = temp;
			}

			// Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
			int c = (i * (i + 1))/2 + j;

			if(c < 0 || c > 54) {
				return NULL;
			}

			// Return a pointer to this memory location
			return &(M_alpha[c]);
		}

		void build(mesh_node *n[10])
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
			scalar psi[10];

			SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);

			zero();

			scalar J_inv[9];

			for(int i = 0; i < NUM_TET_GAUSS_QUAD_POINTS; i++) {
				scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
				SecondOrderFunctions::calc_psi(psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2]);
//				if(i == 0) printf("MassMatrix: det_J = %e\n", det_J);
				add_psi_dot_products(psi, det_J, gauss_points[i].W);
			}


/*
			printf("MASSMATRIXQUADRATIC:\n");
			scalar M[10][10];
			int c = 0;
			for(int j = 0; j < 10; j++) {
				for(int i = 0; i <= j; i++) {
					M[i][j] = M_alpha[c];
					M[j][i] = M_alpha[c];
					c++;
				}
			}

			for(int i = 0; i < 10; i++) {
				for(int j = 0; j < 10; j++) {
					printf("%+e ", M[i][j]);
				}
				printf("\n");
			}
			printf("\n");
*/
		}

		scalar get_M_alpha_value(int i, int j)
		{
			// Poisson matrix is symmetric, so convert any request for an upper triangular element into its
			// corresponding (equivalent) lower triangular element
			if(i < j) {
				int temp = i;
				i = j;
				j = temp;
			}

			// Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
			int c = (i * (i + 1))/2 + j;

			if(c < 0 || c > 54) {
				return -897000;
			}

			// Return a pointer to this memory location
			return M_alpha[c];
		}

	private:
		scalar M_alpha[NUM_ELEMENTS_LOWER_TRIANGULAR_10X10];

                struct tetrahedron_gauss_point
                {
                        scalar W;
                        scalar eta[4];
                };

		void add_psi_dot_products(scalar psi[10], scalar det_J, scalar weight)
		{
			int c = 0;
			for(int i = 0; i < 10; i++) {
				for(int j = 0; j <= i; j++) {
					M_alpha[c] += det_J * weight * psi[i] * psi[j];
					c++;
				}
			}
		}

		void zero()
		{
			for(int i = 0; i < NUM_ELEMENTS_LOWER_TRIANGULAR_10X10; i++) {
				M_alpha[i] = 0;
			}
		}
};


#endif
