#include "NoMassCGSolver.h"

NoMassCGSolver::NoMassCGSolver() {
    num_rows = 0;
    num_nodes = 0;
    epsilon2 = 0;
    i_max = 0;
    preconditioner = NULL;
    r = NULL;
    p = NULL;
    z = NULL;
    q = NULL;
    f = NULL;
}

/* */
NoMassCGSolver::~NoMassCGSolver() {
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] q;
    delete[] f;
    delete[] preconditioner;
    r = NULL;
    p = NULL;
    z = NULL;
    q = NULL;
    f = NULL;
    preconditioner = NULL;
    num_rows = 0;
    num_nodes = 0;
    epsilon2 = 0;
    i_max = 0;
}

/* */
int NoMassCGSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list) {

    this->num_rows = 3 * num_nodes;
    this->num_nodes = num_nodes;
    this->epsilon2 = params->epsilon2;
    this->i_max = params->max_iterations_cg;
    this->one = 1;
    printf("\t\t\tCalculating Sparsity Pattern for a 1st Order Viscosity Matrix\n");
    SparsityPattern sparsity_pattern_viscosity_matrix;
    sparsity_pattern_viscosity_matrix.init(num_rows);

    scalar *mem_loc;
    int n, ni, nj, ni_index, nj_index, i, j;
    matrix3 J;
    for (n = 0; n < num_elements; n++) {
        elem[n].calculate_jacobian(J);
        elem[n].calc_shape_function_derivatives_and_volume(J);
        elem[n].create_viscosity_matrix();
        for (ni = 0; ni < 10; ++ni) {
            for (nj = 0; nj < 10; ++nj) {
                ni_index = 3 * elem[n].n[ni]->index;
                nj_index = 3 * elem[n].n[nj]->index;
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        if (ni < 4 && nj < 4) {
                            mem_loc = &elem[n].viscosity_matrix[ni + 4 * i][nj + 4 * j];
                            sparsity_pattern_viscosity_matrix.register_contribution(ni_index + i, nj_index + j, mem_loc);
                        } else {
                            if (ni == nj && i == j) {
                                if (sparsity_pattern_viscosity_matrix.check_for_contribution(ni_index + i, nj_index + j) == false) {
                                    mem_loc = &one;
                                    sparsity_pattern_viscosity_matrix.register_contribution(ni_index + i, nj_index + j, mem_loc);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (params->do_stokes == 1) {
        for (ni = 0; ni < num_nodes; ++ni) {
            for (nj = 0; nj < 3; ++nj) {
                sparsity_pattern_viscosity_matrix.register_contribution(3 * ni + nj, 3 * ni + nj, &node[ni].stokes_drag);
            }
        }
    }

    // Creating fixed pattern viscosity matrix, but not entering values! Just making the pattern for now
    printf("\t\t\tBuilding Sparsity Pattern for Viscosity Matrix\n");
    V = sparsity_pattern_viscosity_matrix.create_sparse_matrix();

    // create a preconditioner for solving in less iterations
    // Create the jacobi preconditioner matrix (diagonal)
    preconditioner = new scalar[num_rows];

    // create the work vectors necessary for use by the conjugate gradient solver in 'solve'
    r = new vector3[num_nodes];
    p = new vector3[num_nodes];
    z = new vector3[num_nodes];
    q = new vector3[num_nodes];
    f = new vector3[num_nodes];

    return FFEA_OK;
}

/*  */
int NoMassCGSolver::solve(vector3* x) {
    // Complete the sparse viscosity matrix
    V->build();
    V->calc_inverse_diagonal(preconditioner);
    //V->print_dense_to_file(x);
    int i = 0;
    scalar delta_new, delta_old, pTq, alpha;
    delta_new = conjugate_gradient_residual_assume_x_zero(x);

    for (i = 0; i < i_max; i++) {
        pTq = get_alpha_denominator();
        alpha = delta_new / pTq;

        // Update solution and residual
        vec3_add_to_scaled(x, p, alpha, num_nodes);
        vec3_add_to_scaled(r, q, -alpha, num_nodes);

        // Once convergence is achieved, return
        if (residual2() < epsilon2) {
            return FFEA_OK;
        }

        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();
        vec3_scale_and_add(p, z, (delta_new / delta_old), num_nodes);
    }
    // If desired convergence was not reached in the set number of iterations...
    FFEA_ERROR_MESSG("Conjugate gradient solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", i_max);
}

/* */
void NoMassCGSolver::print_matrices(vector3* force) {
    V->print_dense_to_file(force);
}

/* */
scalar NoMassCGSolver::conjugate_gradient_residual_assume_x_zero(vector3 *b) {
    int i = 0;
    scalar delta_new = 0;

#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(b, stderr) reduction(+:delta_new)
#endif
    for (i = 0; i < num_nodes; i++) {
        r[i].x = b[i].x;
        r[i].y = b[i].y;
        r[i].z = b[i].z;
        f[i].x = b[i].x;
        f[i].y = b[i].y;
        f[i].z = b[i].z;
        b[i].x = 0;
        b[i].y = 0;
        b[i].z = 0;
        z[i].x = preconditioner[(3 * i)] * r[i].x;
        z[i].y = preconditioner[(3 * i) + 1] * r[i].y;
        z[i].z = preconditioner[(3 * i) + 2] * r[i].z;
        p[i].x = z[i].x;
        p[i].y = z[i].y;
        p[i].z = z[i].z;
        delta_new += r[i].x * z[i].x + r[i].y * z[i].y + r[i].z * z[i].z;
    }
    return delta_new;
}

/* */
scalar NoMassCGSolver::residual2() {
    int i;
    scalar r2 = 0, f2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(stderr) reduction(+:r2, f2)
#endif
    for (i = 0; i < num_nodes; i++) {
        r2 += r[i].x * r[i].x + r[i].y * r[i].y + r[i].z * r[i].z;
        f2 += f[i].x * f[i].x + f[i].y * f[i].y + f[i].z * f[i].z;
    }
    if (f2 == 0.0) {
        return 0.0;
    } else {
        return r2 / f2;
    }
}

/* */
scalar NoMassCGSolver::modx(vector3 *x) {
    int i;
    scalar r2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(x) reduction(+:r2)
#endif
    for (i = 0; i < num_nodes; i++) {
        r2 += x[i].x * x[i].x + x[i].y * x[i].y + x[i].z * x[i].z;
    }
    return r2;
}

scalar NoMassCGSolver::get_alpha_denominator() {
    // A * p
    V->apply(p, q);
    int i;
    scalar pTq = 0;

    // p^T * A * p
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) reduction(+:pTq)
#endif
    for (i = 0; i < num_nodes; ++i) {
        pTq += p[i].x * q[i].x + p[i].y * q[i].y + p[i].z * q[i].z;
    }

    return pTq;
}

/* */
scalar NoMassCGSolver::parallel_apply_preconditioner() {
    int i;
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) reduction(+:delta_new)
#endif
    for (i = 0; i < num_nodes; i++) {
        z[i].x = preconditioner[(3 * i)] * r[i].x;
        z[i].y = preconditioner[(3 * i) + 1] * r[i].y;
        z[i].z = preconditioner[(3 * i) + 2] * r[i].z;
        delta_new += r[i].x * z[i].x + r[i].y * z[i].y + r[i].z * z[i].z;
    }

    return delta_new;
}

/* */
void NoMassCGSolver::check(vector3 *x) {
    FILE *fout2;
    fout2 = fopen("/localhome/py09bh/output/nomass/cube_viscosity_no_mass.csv", "a");
    int i;
    double temp = 0, temp2 = 0;
    vector3 *temp_vec = new vector3[num_nodes];
    V->apply(x, temp_vec);
    for (i = 0; i < num_nodes; ++i) {
        temp += x[i].x * temp_vec[i].x + x[i].y * temp_vec[i].y + x[i].z * temp_vec[i].z;
        temp2 += x[i].x * f[i].x + x[i].y * f[i].y + x[i].z * f[i].z;
    }
    fprintf(fout2, "%e,%e\n", temp2, fabs(temp - temp2));
    fclose(fout2);

}

/* */



