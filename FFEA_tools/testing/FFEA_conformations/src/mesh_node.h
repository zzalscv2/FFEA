/*
 *	mesh_node.hpp
 *
 */

#ifndef MESH_NODE_H_INCLUDED
#define MESH_NODE_H_INCLUDED

#include <stddef.h>
#include "mat_vec_fns.h"

/*
 * Structure for a mesh_node: the points FEM meshes are built from.
 */
class mesh_node
{
	public:

		mesh_node()
		{
			num_element_contributors = 0;
			force_contributions = NULL;
			pos.x = 0; pos.y = 0; pos.z = 0;
			vel.x = 0; vel.y = 0; vel.z = 0;
			phi = 0;
			index = 0;
			pos_0.x = 0; pos_0.y = 0; pos_0.z = 0;
			stokes_radius = 0;
			stokes_drag = 0;
		}

		~mesh_node()
		{
			delete[] force_contributions;
			pos.x = 0; pos.y = 0; pos.z = 0;
			vel.x = 0; vel.y = 0; vel.z = 0;
			phi = 0;
			index = 0;
			pos_0.x = 0; pos_0.y = 0; pos_0.z = 0;
			num_element_contributors = 0;
			stokes_radius = 0;
			stokes_drag = 0;
		}

		void print()
		{
			printf("pos: %e %e %e\n", pos.x, pos.y, pos.z);
			printf("vel: %e %e %e\n", vel.x, vel.y, vel.z);
		}

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
};

#endif
