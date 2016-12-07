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
