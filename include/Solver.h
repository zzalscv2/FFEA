#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

#include "tetra_element_linear.h"
#include <set>

class Solver {
public:

    /** Provide an empty constructor (to avoid linker problems) */
    Solver() {
    };

    /** Make the destructor virtual (so that the destructor of derived classes will be called too) */
    virtual ~Solver() {
    };

    /**
     * Initialises the solver (by building whatever representation of the mass matrix it needs)
     * using the given node-element connectivity.
     */
    virtual int init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_nodes_list) = 0;

    /**
     * Solves the linear system Mx = f where f is the force vector (should be 'x' on input), and the mass matrix
     * M has already been constructed (in whatever representation) by the init() function for
     * a particular Blob. The solution is written to x.
     */
    virtual int solve(vector3 *x) = 0;

    virtual void apply_matrix(scalar *in, scalar *result) = 0;
};
#endif
