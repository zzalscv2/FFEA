/*
 *	mesh_node.cpp
 *
 */

#include "mesh_node.h"

/*
 * Structure for a mesh_node: the points FEM meshes are built from.
 */
mesh_node::mesh_node() {
    num_element_contributors = 0;
    force_contributions = NULL;
    pos.x = 0;
    pos.y = 0;
    pos.z = 0;
    vel.x = 0;
    vel.y = 0;
    vel.z = 0;
    phi = 0;
    index = 0;
    pos_0.x = 0;
    pos_0.y = 0;
    pos_0.z = 0;
    stokes_radius = 0;
    stokes_drag = 0;
    linear = false;
}

mesh_node::~mesh_node() {
    delete[] force_contributions;
    pos.x = 0;
    pos.y = 0;
    pos.z = 0;
    vel.x = 0;
    vel.y = 0;
    vel.z = 0;
    phi = 0;
    index = 0;
    pos_0.x = 0;
    pos_0.y = 0;
    pos_0.z = 0;
    num_element_contributors = 0;
    stokes_radius = 0;
    stokes_drag = 0;
    linear = false;
}

void mesh_node::move(int direction, scalar dx) {
	switch(direction) {
		case(0):
			pos.x += dx;
			break;
		case(1):
			pos.y += dx;
			break;
		case(2):
			pos.z += dx;
			break;
	}
}

void mesh_node::set_pos(scalar x, scalar y, scalar z) {
	pos.x = x;
	pos.y = y;
	pos.z = z;
}

void mesh_node::print() {
    printf("pos: %e %e %e\n", pos.x, pos.y, pos.z);
    printf("vel: %e %e %e\n", vel.x, vel.y, vel.z);
}

void mesh_node::set_linear() {
    linear = true;
}

bool mesh_node::am_I_linear() {
    return linear;
}