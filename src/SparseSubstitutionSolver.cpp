#include "SparseSubstitutionSolver.h"

SparseSubstitutionSolver::SparseSubstitutionSolver() {
    // Initialise everything to zero
    num_rows = 0;
    L_key = NULL;
    U_key = NULL;
    inverse_diag = NULL;
    L = NULL;
    U = NULL;
}

SparseSubstitutionSolver::~SparseSubstitutionSolver() {
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

int SparseSubstitutionSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list) {

    // Mass matrix will have as many rows as there are nodes in the mesh
    num_rows = num_nodes;

    printf("Allocating memory for mass matrix...\n");
    scalar *mass = new scalar[num_rows * num_rows];
    scalar *mass_LU = new scalar[num_rows * num_rows];
    printf("...done.\n");

    printf("Zeroing...\n");
    for (int i = 0; i < num_rows * num_rows; i++) {
        mass[i] = 0;
        mass_LU[i] = 0;
    }
    printf("...done\n");

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    int is_pinned[num_nodes];
    for (int i = 0; i < num_nodes; i++) {
        is_pinned[i] = 0;
    }
    for (int i = 0; i < num_pinned_nodes; i++) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }

    // build the matrix
    printf("Building the mass matrix...\n");
    int ni, nj;
    for (int n = 0; n < num_elements; n++) {
        // add mass matrix for this element
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (i < 4 && j < 4) {
                    ni = elem[n].n[i]->index;
                    nj = elem[n].n[j]->index;
                    if (is_pinned[ni] == 0 && is_pinned[nj] == 0) {
                        if (i == j) {
                            mass[ni * num_rows + nj] += .1 * elem[n].rho * elem[n].vol_0;
                        } else {
                            mass[ni * num_rows + nj] += .05 * elem[n].rho * elem[n].vol_0;
                        }
                    } else {
                        if (i == j) {
                            mass[ni * num_rows + nj] = 1;
                        }
                    }
                } else {
                    if (i == j) {
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
    for (int k = 0; k < num_rows; k++) {
        for (int i = 0; i < k + 1; i++) {
            scalar sum_ij_kj = 0;
            for (int j = 0; j < i; j++) {
                sum_ij_kj += mass_LU[INDEX(i, j)] * mass_LU[INDEX(k, j)];
            }
            mass_LU[INDEX(k, i)] = (mass[INDEX(k, i)] - sum_ij_kj) / mass_LU[INDEX(i, i)];
        }

        scalar sum_kj_2 = 0;
        for (int j = 0; j < k; j++) {
            sum_kj_2 += mass_LU[INDEX(k, j)] * mass_LU[INDEX(k, j)];
        }
        mass_LU[INDEX(k, k)] = sqrt(mass[INDEX(k, k)] - sum_kj_2);
    }
    delete[] mass;
    printf("...done.\n");

    // Copy the lower matrix into upper
    for (int i = 1; i < num_rows; i++) {
        for (int j = 0; j < i; j++) {
            mass_LU[j * num_rows + i] = mass_LU[i * num_rows + j];
        }
    }

    // Allocate and fill the 'inverse_diag' array
    inverse_diag = new scalar[num_rows];

    for (int i = 0; i < num_rows; i++)
        inverse_diag[i] = 1.0 / mass_LU[i * num_rows + i];

    // Allocate the 'key' arrays
    L_key = new int[num_rows];
    U_key = new int[num_rows];

    // Build the lower triangular matrix key
    int total_L = 0;
    for (int i = 0; i < num_rows; i++) {
        for (int j = num_rows - 1; j >= i; j--) {
            if (mass_LU[i * num_rows + j] != 0) {
                total_L += (j - i);
                L_key[i] = (j - i);
                break;
            }
        }
    }

    // Build the upper triangular matrix key
    total_entries_in_U = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = 0; j <= i; j++)
            if (mass_LU[i * num_rows + j] != 0) {
                total_entries_in_U += (i - j);
                U_key[i] = (i - j);
                break;
            }

    // Allocate the off-diagonal entry arrays
    L = new scalar[total_L];
    U = new scalar[total_entries_in_U];

    // Fill up the L and U sparse triangular matrices
    int off_diag_data_index = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = i + 1; j <= i + L_key[i]; j++) {
            L[off_diag_data_index] = mass_LU[i * num_rows + j];
            off_diag_data_index++;
        }
    off_diag_data_index = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = i - U_key[i]; j < i; j++) {
            U[off_diag_data_index] = mass_LU[i * num_rows + j];
            off_diag_data_index++;
        }

    delete[] mass_LU;

    return FFEA_OK;
}

int SparseSubstitutionSolver::solve(vector3 *x) {
    int i, j, index;

    // Forward substitution step Ly = b :
    index = 0;
    for (i = 0; i < num_rows; i++) {

        x[i].x *= inverse_diag[i];
        x[i].y *= inverse_diag[i];
        x[i].z *= inverse_diag[i];

        for (j = 0; j < L_key[i]; j++) {
            x[i + j + 1].x -= x[i].x * L[j + index];
            x[i + j + 1].y -= x[i].y * L[j + index];
            x[i + j + 1].z -= x[i].z * L[j + index];
        }

        index += L_key[i];
    }

    // Backward substitution step Ux = y :
    index = total_entries_in_U - 1;
    for (i = num_rows - 1; i >= 0; i--) {

        x[i].x *= inverse_diag[i];
        x[i].y *= inverse_diag[i];
        x[i].z *= inverse_diag[i];

        for (j = 0; j < U_key[i]; j++) {
            x[i - j - 1].x -= x[i].x * U[index - j];
            x[i - j - 1].y -= x[i].y * U[index - j];
            x[i - j - 1].z -= x[i].z * U[index - j];
        }

        index -= U_key[i];
    }

    return FFEA_OK;
}

void SparseSubstitutionSolver::apply_matrix(scalar *in, scalar *result) {
}

