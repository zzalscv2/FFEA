#ifndef CG_SOLVER_H_INCLUDED
#define CG_SOLVER_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "SparseMatrixFixedPattern.h"

class CG_solver
{
	public:
		CG_solver();
		~CG_solver();

		int init(int N, scalar tol, int max_num_iterations);

		int solve(SparseMatrixFixedPattern *A, scalar *x, scalar *b);

		int solve(SparseMatrixFixedPattern *A, scalar *x, scalar *b, int num_iterations);

	private:

		// The number of unknowns
		int N;

		// The convergence tolerance threshold
		scalar tol;

		// Maximum number of iterations before giving up
		int max_num_iterations;

		// The preconditioner matrix (inverse of the diagonal)
		scalar *inv_M;

		// Vectors needed for use by conjugate gradient solver
		scalar *d, *r, *q, *s;

		scalar conjugate_gradient_residual(SparseMatrixFixedPattern *A, scalar *x, scalar *b);

		scalar residual2();

		void parallel_vector_add_self(scalar *v1, scalar a, scalar *v2);

		void parallel_vector_add(scalar *v1, scalar a, scalar *v2);

		scalar parallel_apply_preconditioner();

		void zero(scalar *v);

		/* Returns the dot product of vectors a and b, of length N */
		scalar dot(scalar *a, scalar *b);
};

#endif
