#ifndef BICGSTAB_SOLVER_H_INCLUDED
#define BICGSTAB_SOLVER_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "SparseMatrixUnknownPattern.h"

class BiCGSTAB_solver {
public:
    BiCGSTAB_solver();

    ~BiCGSTAB_solver();

    int init(int N, scalar tol, int max_num_iterations);

    int solve(SparseMatrixUnknownPattern *A, scalar *x, scalar *b);

    int solve(SparseMatrixUnknownPattern *A, scalar *x, scalar *b, int num_iterations);

private:

    /// The number of unknowns
    int N;

    /// The convergence tolerance threshold
    scalar tol;

    /// Maximum number of iterations before giving up
    int max_num_iterations;

    /// The inverse of the preconditioner matrix
    scalar *inv_M;

    //@{
    /// The residual vectors
    scalar *r, *r_hat;
    //@}

    //@{
    /// Other necessary vectors
    scalar *p, *p_hat, *q, *s, *s_hat, *t;
    //@}

    /** Calculates the residual vector, r, for the matrix equation Ax = b. Specifically, r = b - Ax. */
    void get_residual_vector(scalar *r, scalar *b, SparseMatrixUnknownPattern *A, scalar *x, int N);

    /** Copies the contents of vector b into vector a (a <- b) */
    void copy_vector(scalar *a, scalar *b, int N);

    /** Returns the dot product of vectors a and b, of length N */
    scalar dot(scalar *a, scalar *b, int N);

    /** Calculates y = Mx for diagonal matrix and vectors of dimension N */
    void apply_diagonal_matrix(scalar *y, scalar *M, scalar *x, int N);

    /** Sets the given vector (length N) to zero */
    void zero(scalar *x, int N);

    /** Carries out the operation x = y + c*z, where c is a scalar, and x, y and z are vectors of length N. */
    void scalar_vector_add(scalar *x, scalar *y, scalar c, scalar *z, int N);

    /** Carries out the operation w = x + a * (y + b * z) */
    // p = r + beta(p - omega * v)
    void complicated_machine(scalar *w, scalar *x, scalar a, scalar *y, scalar b, scalar *z, int N);
};

#endif
