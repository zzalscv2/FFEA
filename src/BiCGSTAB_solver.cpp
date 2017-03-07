// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

#include "BiCGSTAB_solver.h"

BiCGSTAB_solver::BiCGSTAB_solver() {
    N = 0;
    tol = 0;
    max_num_iterations = 0;

    inv_M = NULL;
    r = NULL;
    r_hat = NULL;
    p = NULL;
    p_hat = NULL;
    q = NULL;
    s = NULL;
    s_hat = NULL;
    t = NULL;
}

BiCGSTAB_solver::~BiCGSTAB_solver() {
    delete[] inv_M;
    delete[] r;
    delete[] r_hat;
    delete[] p;
    delete[] p_hat;
    delete[] q;
    delete[] s;
    delete[] s_hat;
    delete[] t;

    inv_M = NULL;
    r = NULL;
    r_hat = NULL;
    p = NULL;
    p_hat = NULL;
    q = NULL;
    s = NULL;
    s_hat = NULL;
    t = NULL;

    N = 0;
    tol = 0;
    max_num_iterations = 0;
}

int BiCGSTAB_solver::init(int N, scalar tol, int max_num_iterations) {
    this->N = N;
    this->tol = tol;
    this->max_num_iterations = max_num_iterations;

    // Allocate all required memory
    inv_M = new(std::nothrow) scalar[N];
    r = new(std::nothrow) scalar[N];
    r_hat = new(std::nothrow) scalar[N];
    p = new(std::nothrow) scalar[N];
    p_hat = new(std::nothrow) scalar[N];
    q = new(std::nothrow) scalar[N];
    s = new(std::nothrow) scalar[N];
    s_hat = new(std::nothrow) scalar[N];
    t = new(std::nothrow) scalar[N];

    // Check that memory has been allocated
    if (inv_M == NULL ||
            r == NULL ||
            r_hat == NULL ||
            p == NULL ||
            p_hat == NULL ||
            q == NULL ||
            s == NULL ||
            s_hat == NULL ||
            t == NULL) {
        FFEA_ERROR_MESSG("While initialising BiCGSTAB_solver, could not allocate memory for vectors.\n");
    }

    // Initialise all vectors to zero
    zero(inv_M, N);
    zero(r, N);
    zero(r_hat, N);
    zero(p, N);
    zero(p_hat, N);
    zero(q, N);
    zero(s, N);
    zero(s_hat, N);
    zero(t, N);

    return FFEA_OK;
}

int BiCGSTAB_solver::solve(SparseMatrixUnknownPattern *A, scalar *x, scalar *b) {
    scalar rho_last = 1, alpha = 1, omega = 1, rho, beta;

    // Get the inverse of the diagonal of the matrix to use as a preconditioner
    A->calc_inverse_diagonal(inv_M);

    // Get the initial residual vector using a first guess at the solution x (r_0 = b - A x_0)
    get_residual_vector(r, b, A, x, N);

    // r_hat_0 <- r_0
    copy_vector(r_hat, r, N);

    // initialise q vector to zero
    zero(q, N);

    // loop till convergence (or a threshold maximum number of iterations have elapsed -> convergence failure)
    int i;
    for (i = 0; i < max_num_iterations; i++) {

        // rho_i = dot(r_hat_0, r_i_minus_1);
        rho = dot(r_hat, r, N);

        if (rho == 0) {
            FFEA_ERROR_MESSG("In BiCGSTAB_solver solve(), rho is zero. Solver stopping on iteration %d.\n", i);
        }

        // beta = (rho_i/rho_i_minus_1) (alpha/omega_i_minus_1)
        beta = (rho / rho_last) * (alpha / omega);

        // p = r + beta * (p - omega * q);
        complicated_machine(p, r, beta, p, -omega, q, N);

        // apply preconditioner to p (p_hat = inv_M p_i)
        apply_diagonal_matrix(p_hat, inv_M, p, N);

        // q = A p_hat
        A->apply(p_hat, q);

        // alpha = rho_i / dot(r_hat_0, q)
        alpha = rho / dot(r_hat, q, N);

        // s = r - alpha * q;
        scalar_vector_add(s, r, -alpha, q, N);

        // x = x + alpha * p_hat
        scalar_vector_add(x, x, alpha, p_hat, N);

        // Check for convergence
        if (dot(s, s, N) < tol) {
            //			printf("Convergence reached on iteration %d.\n", i);
            return FFEA_OK;
        }

        // apply preconditioner to s (s_hat = inv_M s)
        apply_diagonal_matrix(s_hat, inv_M, s, N);

        // t = A s_hat
        A->apply(s_hat, t);

        // omega_i = dot(t, s) / dot(t, t)
        omega = dot(t, s, N) / dot(t, t, N);

        // x = x + omega * s_hat;
        scalar_vector_add(x, x, omega, s_hat, N);

        // r = s - omega * t;
        scalar_vector_add(r, s, -omega, t, N);

        //		printf("BiCGSTAB_solver: %f\n", dot(s,s,N));
    }

    FFEA_ERROR_MESSG("Bi-Conjugate Gradient Stabilised solver could not converge in max_num_iterations.\n");
}

int BiCGSTAB_solver::solve(SparseMatrixUnknownPattern *A, scalar *x, scalar *b, int num_iterations) {
    scalar rho_last = 1, alpha = 1, omega = 1, rho, beta;

    // Get the inverse of the diagonal of the matrix to use as a preconditioner
    A->calc_inverse_diagonal(inv_M);

    // Get the initial residual vector using a first guess at the solution x (r_0 = b - A x_0)
    get_residual_vector(r, b, A, x, N);

    // r_hat_0 <- r_0
    copy_vector(r_hat, r, N);

    // initialise q vector to zero
    zero(q, N);

    // loop till convergence (or a threshold maximum number of iterations have elapsed -> convergence failure)
    int i;
    for (i = 0; i < num_iterations; i++) {

        // rho_i = dot(r_hat_0, r_i_minus_1);
        rho = dot(r_hat, r, N);

        if (rho == 0) {
            FFEA_ERROR_MESSG("In BiCGSTAB_solver solve(), rho is zero. Solver stopping on iteration %d.\n", i);
        }

        // beta = (rho_i/rho_i_minus_1) (alpha/omega_i_minus_1)
        beta = (rho / rho_last) * (alpha / omega);

        // p = r + beta * (p - omega * q);
        complicated_machine(p, r, beta, p, -omega, q, N);

        // apply preconditioner to p (p_hat = inv_M p_i)
        apply_diagonal_matrix(p_hat, inv_M, p, N);

        // q = A p_hat
        A->apply(p_hat, q);

        // alpha = rho_i / dot(r_hat_0, q)
        alpha = rho / dot(r_hat, q, N);

        // s = r - alpha * q;
        scalar_vector_add(s, r, -alpha, q, N);

        // x = x + alpha * p_hat
        scalar_vector_add(x, x, alpha, p_hat, N);

        // apply preconditioner to s (s_hat = inv_M s)
        apply_diagonal_matrix(s_hat, inv_M, s, N);

        // t = A s_hat
        A->apply(s_hat, t);

        // omega_i = dot(t, s) / dot(t, t)
        omega = dot(t, s, N) / dot(t, t, N);

        // x = x + omega * s_hat;
        scalar_vector_add(x, x, omega, s_hat, N);

        // r = s - omega * t;
        scalar_vector_add(r, s, -omega, t, N);

        //		printf("BiCGSTAB_solver: %f\n", dot(s,s,N));
    }

    return FFEA_OK;
}

void BiCGSTAB_solver::get_residual_vector(scalar *r, scalar *b, SparseMatrixUnknownPattern *A, scalar *x, int N) {
    // Ax
    A->apply(x, r);

    // r = b - Ax
    int i;
    for (i = 0; i < N; i++)
        r[i] = b[i] - r[i];
}

/* Copies the contents of vector b into vector a (a <- b) */
void BiCGSTAB_solver::copy_vector(scalar *a, scalar *b, int N) {
    int i;
    for (i = 0; i < N; i++)
        a[i] = b[i];
}

/* Returns the dot product of vectors a and b, of length N */
scalar BiCGSTAB_solver::dot(scalar *a, scalar *b, int N) {
    scalar result = 0;
    int i;
    for (i = 0; i < N; i++)
        result += a[i] * b[i];

    return result;
}

/* Calculates y = Mx for diagonal matrix and vectors of dimension N */
void BiCGSTAB_solver::apply_diagonal_matrix(scalar *y, scalar *M, scalar *x, int N) {
    int i;
    for (i = 0; i < N; i++)
        y[i] = M[i] * x[i];
}

/* Sets the given vector (length N) to zero */
void BiCGSTAB_solver::zero(scalar *x, int N) {
    int i;
    for (i = 0; i < N; i++)
        x[i] = 0;
}

/* Carries out the operation x = y + c*z, where c is a scalar, and x, y and z are vectors of length N. */
void BiCGSTAB_solver::scalar_vector_add(scalar *x, scalar *y, scalar c, scalar *z, int N) {
    int i;
    for (i = 0; i < N; i++)
        x[i] = y[i] + c * z[i];
}

/* Carries out the operation w = x + a * (y + b * z) */
// p = r + beta(p - omega * v)

void BiCGSTAB_solver::complicated_machine(scalar *w, scalar *x, scalar a, scalar *y, scalar b, scalar *z, int N) {
    int i;
    for (i = 0; i < N; i++)
        w[i] = x[i] + a * (y[i] + b * z[i]);
}
