#ifndef NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED
#define NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED

#include <math.h>
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "LinkedListCube.h"
#include "Face.h"

class NearestNeighbourLinkedListCube : public LinkedListCube<Face>
{
	public:

		/* Build the nearest neighbour look up cube given the spatial cell size */
		int build_nearest_neighbour_lookup(scalar h)
		{
			// Clear the grid
			clear();

			// Loop through each Face in the pool, calculate which cell of the
			// grid it belongs in based on its centroid position, and add it to
			// the linked list stack on that cell
			int i;
			for(i = 0; i < num_nodes_in_pool; i++) {

				// calculate which cell the face belongs in
				int x = (int)floor(pool[i].obj->centroid.x/h);
				int y = (int)floor(pool[i].obj->centroid.y/h);
				int z = (int)floor(pool[i].obj->centroid.z/h);
/*
				// If face centroid is out of bounds of the box, add its node to the nearest cell on the edge
				// of the box. This is necessary in the case of PBC since an entire blob is moved only if its centre
				// of mass exceeds the bounds of the box
				if(x < 0) {
					x = 0;
				}
				if(x > N - 1) {
					x -= N - 1;
				}
				if(y < 0) {
					y = 0;
				}
				if(y > N - 1) {
					y = N - 1;
				}
				if(z < 0) {
					z = 0;
				}
				if(z > N - 1) {
					z = N - 1;
				}
*/
				// attempt to add the node to the cell
				if(add_node_to_stack(i, x, y, z) == FFEA_ERROR) {
					FFEA_ERROR_MESSG("Error when trying to add node %d to nearest neighbour stack at (%d %d %d)\n", i, x, y, z);
				}
			}

			return FFEA_OK;
		}
};

#endif
