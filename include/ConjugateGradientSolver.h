#ifndef CONJUGATEGRADIENTSOLVER_HPP_INCLUDED
#define CONJUGATEGRADIENTSOLVER_HPP_INCLUDED


#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "ConjugateGradientSolver.h"
#include "SparseMatrixTypes.h"

class ConjugateGradientSolver : public Solver {
public:
    /** Constructor */
    ConjugateGradientSolver();

    /** Destructor */
    ~ConjugateGradientSolver();

    /** Builds the sparse mass matrix and allocates the various work vectors required for conjugate gradient */
    int init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list);

    /** Applies conjugate gradient with a Jacobi preconditioner to solve the system Mx = f */
    int solve(vector3 *x);

    /** Applies the mass matrix to the given vector, 'in', putting the result in 'result'*/
    void apply_matrix(scalar *in, scalar *result);

private:

    /** Error tolerance threshold (normalised and squared) to determine when solution has converged */
    scalar epsilon2;

    /** Maximum number of iterations the solver should use before giving up (as solution is not converging) */
    int i_max;

    /** Number of rows in original matrix */
    int num_rows;

    /** Array of all non-zero entries comprising the mass matrix, in the order they appear in the matrix
     * when scanned left to right across rows first (and columns secondary)
     */
    sparse_entry *entry;

    /**
     * Array of size (num_nodes+1) containing the index (in the sparse_entry array above)
     * of the first non-zero entry in each row of the original matrix.
     * The final entry in this vector, key[num_nodes], is the total number of non-zero entries
     * in the entire matrix.
     */
    int *key;

    /** Jacobi preconditioner (inverse of the mass matrix diagonal) */
    scalar *preconditioner;

    /** Work vectors */
    vector3 *d, *r, *q, *s, *f;

    /* */
    scalar conjugate_gradient_residual_assume_x_zero(vector3 *b);

    /* */
    scalar parallel_sparse_matrix_apply();

    void parallel_vector_add_self(vector3 *v1, scalar a, vector3 *v2, int vec_size);

    void parallel_vector_add(vector3 *v1, scalar a, vector3 *v2, int vec_size);

    /* */
    scalar parallel_apply_preconditioner();

    /* */
    scalar residual2();
};

#endif
