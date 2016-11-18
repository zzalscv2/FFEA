#include "MassLumpedSolver.h"

/**/
MassLumpedSolver::MassLumpedSolver() {
    num_rows = 0;
    inv_M = NULL;
}

/* */
MassLumpedSolver::~MassLumpedSolver() {
    delete[] inv_M;
    inv_M = NULL;
    num_rows = 0;
}

/* */
int MassLumpedSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list) {
    int n, i, ni;

    // Store the number of rows, error threshold (stopping criterion for solver) and max
    // number of iterations, on this Solver (these quantities will be used a lot)
    this->num_rows = num_nodes;
    inv_M = new scalar[num_rows];

    if (inv_M == NULL) {
        FFEA_ERROR_MESSG("could not allocate inv_M\n");
    }

    for (i = 0; i < num_rows; i++) {
        inv_M[i] = 0;
    }
    for (n = 0; n < num_elements; n++) {
        // add mass contribution for this element
        for (i = 0; i < 10; i++) {
            if (i < 4) {
                ni = elem[n].n[i]->index;
                inv_M[ni] += .25 * elem[n].rho * elem[n].vol_0;
            } else {
                ni = elem[n].n[i]->index;
                inv_M[ni] = 1;
            }
        }
    }

    // set elements corresponding to unmovable 'pinned' nodes to 1
    for (i = 0; i < num_pinned_nodes; i++) {
        inv_M[pinned_nodes_list[i]] = 1.0;
    }

    // inverse
    for (i = 0; i < num_rows; i++) {
        inv_M[i] = 1.0 / inv_M[i];
    }

    return FFEA_OK;
}

/* */
int MassLumpedSolver::solve(vector3 *x) {
    int i = 0;
    for (i = 0; i < num_rows; i++) {
        x[i].x = x[i].x * inv_M[i];
        x[i].y = x[i].y * inv_M[i];
        x[i].z = x[i].z * inv_M[i];
    }
    return FFEA_OK;
}

/* */
void MassLumpedSolver::apply_matrix(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] *= 1.0 / inv_M[i];
    }
}

