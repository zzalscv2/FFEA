#include "CG_solver.h"

		CG_solver::CG_solver()
		{
			N = 0;
			tol = 0;
			max_num_iterations = 0;
			inv_M = NULL;
			d = NULL;
			r = NULL;
			q = NULL;
			s = NULL;
		}

		CG_solver::~CG_solver()
		{
			delete[] inv_M;
			delete[] d;
			delete[] r;
			delete[] q;
			delete[] s;

			N = 0;
			tol = 0;
			max_num_iterations = 0;
			inv_M = NULL;
			d = NULL;
			r = NULL;
			q = NULL;
			s = NULL;
		}

		int CG_solver::init(int N, scalar tol, int max_num_iterations)
		{
			this->N = N;
			this->tol = tol;
			this->max_num_iterations = max_num_iterations;

			// Allocate all required memory
			inv_M = new scalar[N];
			d = new scalar[N];
			r = new scalar[N];
			q = new scalar[N];
			s = new scalar[N];

			// Check that memory has been allocated
			if(	inv_M == NULL	||
				d == NULL	||
				r == NULL	||
				q == NULL	||
				s == NULL )
			{
				FFEA_ERROR_MESSG("While initialising CG_solver, could not allocate memory for vectors.\n");
			}

			// Initialise all vectors to zero
			zero(inv_M);
			zero(d);
			zero(r);
			zero(q);
			zero(s);

			return FFEA_OK;
		}

		int CG_solver::solve(SparseMatrixFixedPattern *A, scalar *x, scalar *b)
		{
			// Get the preconditioner matrix (inverse of the matrix diagonal)
			A->calc_inverse_diagonal(inv_M);

			scalar delta_new, delta_old, alpha;

			// Get the residual vector for Ax = b with current x
			delta_new = conjugate_gradient_residual(A, x, b);

			for(int i = 0; i < max_num_iterations; i++) {

				// Once convergence is achieved, return
				if(residual2() < tol) {
//					printf("CG_solver: Convergence reached on iteration %d\n", i);
					return FFEA_OK;
				}

				// q = A * d
				A->apply(d, q);

				alpha = delta_new/dot(d, q);

				parallel_vector_add_self(x, alpha, d);
				parallel_vector_add_self(r, -alpha, q);

				delta_old = delta_new;
				delta_new = parallel_apply_preconditioner();

				parallel_vector_add(d, (delta_new/delta_old), s);

			}

			// If desired convergence was not reached in the set number of iterations...
			FFEA_ERROR_MESSG("CG_solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", max_num_iterations);
		}

		int CG_solver::solve(SparseMatrixFixedPattern *A, scalar *x, scalar *b, int num_iterations)
		{
			// Get the preconditioner matrix (inverse of the matrix diagonal)
			A->calc_inverse_diagonal(inv_M);

			scalar delta_new, delta_old, alpha;

			// Get the residual vector for Ax = b with current x
			delta_new = conjugate_gradient_residual(A, x, b);

			for(int i = 0; i < num_iterations; i++) {

				// q = A * d
				A->apply(d, q);

				alpha = delta_new/dot(d, q);

				parallel_vector_add_self(x, alpha, d);
				parallel_vector_add_self(r, -alpha, q);

				delta_old = delta_new;
				delta_new = parallel_apply_preconditioner();

				parallel_vector_add(d, (delta_new/delta_old), s);

			}

			return FFEA_OK;
		}

		scalar CG_solver::conjugate_gradient_residual(SparseMatrixFixedPattern *A, scalar *x, scalar *b)
		{
			// Ax
			A->apply(x, r);

			// r = b - Ax
			// d = inv_M * r
			// delta_new = r . d
			scalar delta_new = 0;
			for(int i = 0; i < N; i++) {
				r[i] = b[i] - r[i];
				d[i] = inv_M[i] * r[i];
				delta_new += r[i] * d[i];
			}

			return delta_new;
		}

		scalar CG_solver::residual2()
		{
			int i;
			scalar r2 = 0;
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) private(i) reduction(+:r2)
			#endif
			for(i = 0; i < N; i++) {
				r2 += r[i] * r[i];
			}

			return r2;
		}

		void CG_solver::parallel_vector_add_self(scalar *v1, scalar a, scalar *v2)
		{
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) shared(v1, a, v2)
			#endif
			for(int i = 0; i < N; i++) {
				v1[i] += a * v2[i];
			}
		}

		void CG_solver::parallel_vector_add(scalar *v1, scalar a, scalar *v2)
		{
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) shared(v1, a, v2)
			#endif
			for(int i = 0; i < N; i++) {
				v1[i] = v2[i] + a * v1[i];
			}
		}

		scalar CG_solver::parallel_apply_preconditioner()
		{
			scalar delta_new = 0;
			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				#pragma omp parallel for default(none) reduction(+:delta_new)
			#endif
			for(int i = 0; i < N; i++) {
				s[i] = inv_M[i] * r[i];
				delta_new += r[i] * s[i];
			}

			return delta_new;
		}

		void CG_solver::zero(scalar *v)
		{
			for(int i = 0; i < N; i++) {
				v[i] = 0;
			}
		}

		scalar CG_solver::dot(scalar *a, scalar *b)
		{
			scalar result = 0;
			for(int i = 0; i < N; i++) {
				result += a[i] * b[i];
			}

			return result;
		}
