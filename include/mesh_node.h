/*
 *	mesh_node.hpp
 *
 */

#ifndef MESH_NODE_H_INCLUDED
#define MESH_NODE_H_INCLUDED

#include <stddef.h>
#include <stdio.h>
#include "mat_vec_fns.h"

/*
 * Structure for a mesh_node: the points FEM meshes are built from.
 */
class mesh_node {
public:

    mesh_node();
    ~mesh_node();

    void move(int direction, scalar dx);

    void set_pos(scalar x, scalar y, scalar z);

    void print();

    /* Position of node */
    vector3 pos;

    /* Velocity of node */
    vector3 vel;

    /* Electrostatic potential at this node */
    scalar phi;

    int num_element_contributors;

    /* An array of pointers to contributions to the total force on this node. There should be one
     * contribution from each element this node is a part of (so the length will be num_element_contributors).
     */
    vector3 **force_contributions;

    /* Required for some general matrix constructions in which we need to know this node's 'index' in the node vector */
    int index;

    /* Equilibrium position of nodes (for RMSD calculations) */
    vector3 pos_0;

    /* Charge density on this node */
    scalar rho;

    /* Stokes radius of this node */
    scalar stokes_radius;

    /* The drag due to stokes on this node, not including velocity */
    scalar stokes_drag;

    /* Stores whether or not this node is linear (as the order is surface - interior, not linear - secondary) */
    bool linear;

    void set_linear();
    bool am_I_linear();
};

#endif
