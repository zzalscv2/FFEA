#include "ConjugateGradientSolver.h"

ConjugateGradientSolver::ConjugateGradientSolver() {
    num_rows = 0;
    epsilon2 = 0;
    i_max = 0;
    key = NULL;
    entry = NULL;
    preconditioner = NULL;
    d = NULL;
    r = NULL;
    q = NULL;
    s = NULL;
    f = NULL;
}

ConjugateGradientSolver::~ConjugateGradientSolver() {
    delete[] key;
    delete[] entry;
    delete[] preconditioner;
    delete[] d;
    delete[] r;
    delete[] q;
    delete[] s;
    delete[] f;
    key = NULL;
    entry = NULL;
    preconditioner = NULL;
    d = NULL;
    r = NULL;
    q = NULL;
    s = NULL;
    f = NULL;
    num_rows = 0;
    epsilon2 = 0;
    i_max = 0;
}

int ConjugateGradientSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list) {
    int n, i, j, ni, nj;

    // Store the number of rows, error threshold (stopping criterion for solver) and max
    // number of iterations, on this Solver (these quantities will be used a lot)
    this->num_rows = num_nodes;
    this->epsilon2 = params->epsilon2;
    this->i_max = params->max_iterations_cg;

    printf("\t\tAttempting to allocate %d scalars for mass_LU...\n", num_rows * num_rows);
    scalar *mass_LU = new scalar[num_rows * num_rows];
    printf("\t\t...success.\n");

    if (mass_LU == NULL) {
        FFEA_ERROR_MESSG("Could not allocate mass_LU\n");
    }

    printf("\t\tZeroing...\n");
    for (i = 0; i < num_rows * num_rows; i++) {
        mass_LU[i] = 0;
    }
    printf("\t\t...done\n");

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    int is_pinned[num_nodes];
    for (i = 0; i < num_nodes; i++) {
        is_pinned[i] = 0;
    }
    for (i = 0; i < num_pinned_nodes; i++) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }

    // build the matrix
    printf("\t\tBuilding the mass matrix...\n");
    for (n = 0; n < num_elements; n++) {
        // add mass matrix for this element
        for (i = 0; i < 10; i++) {
            for (j = 0; j < 10; j++) {
                if (i < 4 && j < 4) {
                    ni = elem[n].n[i]->index;
                    nj = elem[n].n[j]->index;
                    if (is_pinned[ni] == 0 && is_pinned[nj] == 0) {
                        if (i == j) {
                            mass_LU[ni * num_rows + nj] += .1 * elem[n].rho * elem[n].vol_0;
                        } else {
                            mass_LU[ni * num_rows + nj] += .05 * elem[n].rho * elem[n].vol_0;
                        }
                    } else {
                        if (i == j) {
                            mass_LU[ni * num_rows + nj] = 1;
                        }
                    }
                } else {
                    if (i == j) {
                        ni = elem[n].n[i]->index;
                        nj = elem[n].n[j]->index;
                        mass_LU[ni * num_rows + nj] = 1;
                    }
                }
            }
        }
    }
    printf("\t\t...done\n");

    // Allocate memory for and initialise 'key' array
    key = new int[num_rows + 1];

    for (i = 0; i < num_rows; i++) {
        key[i] = 0;
    }

    // Get the number of non-zero entries in each row
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_rows; j++) {
            if (mass_LU[i * num_rows + j] != 0) {
                key[i]++;
            }
        }
    }

    // Sum up all the non-zero totals of each row to get the total number of non-zero entries
    // in the whole matrix, constructing the 'key' in the process
    int total_non_zeros = 0, last;
    for (i = 0; i < num_rows; i++) {
        last = key[i];
        key[i] = total_non_zeros;
        total_non_zeros += last;
    }
    key[num_rows] = total_non_zeros;

    // Allocate memory for the 'entry' array
    printf("\t\tNum non zero entries = %d.\n\t\tAllocating space...\n", total_non_zeros);
    entry = new sparse_entry[total_non_zeros];
    printf("\t\t...done.\n");

    // Fill the 'entry' array
    int entry_index = 0;
    scalar val;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_rows; j++) {
            val = mass_LU[i * num_rows + j];
            if (val != 0) {
                entry[entry_index].val = val;
                entry[entry_index].column_index = j;
                entry_index++;
            }
        }
    }

    // Create the jacobi preconditioner matrix (diagonal)
    preconditioner = new scalar[num_rows];
    for (i = 0; i < num_rows; i++)
        preconditioner[i] = 1.0 / mass_LU[i * num_rows + i];

    // create the work vectors necessary for use by the conjugate gradient solver
    d = new vector3[num_rows];
    r = new vector3[num_rows];
    q = new vector3[num_rows];
    s = new vector3[num_rows];
    f = new vector3[num_rows];

    delete[] mass_LU;

    return FFEA_OK;
}

int ConjugateGradientSolver::solve(vector3 *x) {
    int i = 0;
    scalar delta_new, delta_old, dTq, alpha;
    delta_new = conjugate_gradient_residual_assume_x_zero(x);
    for (i = 0; i < i_max; i++) {

        // Once convergence is achieved, return
        if (residual2() < epsilon2) {
            return FFEA_OK;
        }

        dTq = parallel_sparse_matrix_apply();

        alpha = delta_new / dTq;
        parallel_vector_add_self(x, alpha, d, num_rows);
        parallel_vector_add_self(r, -alpha, q, num_rows);

        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();

        parallel_vector_add(d, (delta_new / delta_old), s, num_rows);
    }

    // If desired convergence was not reached in the set number of iterations...
    FFEA_ERROR_MESSG("Conjugate gradient solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", i_max);
}

void ConjugateGradientSolver::apply_matrix(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] = 0;
        for (int j = key[i]; j < key[i + 1]; j++) {
            result[i] += entry[j].val * in[entry[j].column_index];
        }
    }
}

/* */
scalar ConjugateGradientSolver::conjugate_gradient_residual_assume_x_zero(vector3 *b) {
    int i;
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(b) reduction(+:delta_new)
#endif
    for (i = 0; i < num_rows; i++) {
        r[i].x = b[i].x;
        r[i].y = b[i].y;
        r[i].z = b[i].z;
        f[i].x = b[i].x;
        f[i].y = b[i].y;
        f[i].z = b[i].z;
        b[i].x = 0;
        b[i].y = 0;
        b[i].z = 0;
        d[i].x = preconditioner[i] * r[i].x;
        d[i].y = preconditioner[i] * r[i].y;
        d[i].z = preconditioner[i] * r[i].z;
        delta_new += r[i].x * d[i].x + r[i].y * d[i].y + r[i].z * d[i].z;
    }

    return delta_new;
}

/* */
scalar ConjugateGradientSolver::parallel_sparse_matrix_apply() {
    int i, j;
    scalar dTq = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i, j) reduction(+:dTq)
#endif
    for (i = 0; i < num_rows; i++) {
        q[i].x = 0;
        q[i].y = 0;
        q[i].z = 0;
        for (j = key[i]; j < key[i + 1]; j++) {
            q[i].x += entry[j].val * d[entry[j].column_index].x;
            q[i].y += entry[j].val * d[entry[j].column_index].y;
            q[i].z += entry[j].val * d[entry[j].column_index].z;
        }

        dTq += d[i].x * q[i].x + d[i].y * q[i].y + d[i].z * q[i].z;
    }

    return dTq;
}

void ConjugateGradientSolver::parallel_vector_add_self(vector3 *v1, scalar a, vector3 *v2, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, a, v2, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i].x += a * v2[i].x;
        v1[i].y += a * v2[i].y;
        v1[i].z += a * v2[i].z;
    }
}

void ConjugateGradientSolver::parallel_vector_add(vector3 *v1, scalar a, vector3 *v2, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, a, v2, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i].x = v2[i].x + a * v1[i].x;
        v1[i].y = v2[i].y + a * v1[i].y;
        v1[i].z = v2[i].z + a * v1[i].z;
    }
}

/* */
scalar ConjugateGradientSolver::parallel_apply_preconditioner() {
    int i;
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) reduction(+:delta_new)
#endif
    for (i = 0; i < num_rows; i++) {
        s[i].x = preconditioner[i] * r[i].x;
        s[i].y = preconditioner[i] * r[i].y;
        s[i].z = preconditioner[i] * r[i].z;
        delta_new += r[i].x * s[i].x + r[i].y * s[i].y + r[i].z * s[i].z;
    }

    return delta_new;
}

/* */
scalar ConjugateGradientSolver::residual2() {
    int i;
    scalar r2 = 0, f2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(stderr) reduction(+:r2, f2)
#endif
    for (i = 0; i < num_rows; i++) {
        r2 += r[i].x * r[i].x + r[i].y * r[i].y + r[i].z * r[i].z;
        f2 += f[i].x * f[i].x + f[i].y * f[i].y + f[i].z * f[i].z;
    }
    if (f2 == 0.0) {
        return 0.0;
    } else {
        return r2 / f2;
    }
}

