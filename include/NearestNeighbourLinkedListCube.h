#ifndef NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED
#define NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED

#include <math.h>
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "LinkedListCube.h"
#include "Face.h"

class NearestNeighbourLinkedListCube : public LinkedListCube<Face> {
public:
    /* Build the nearest neighbour look up cube given the spatial cell size */
    int build_nearest_neighbour_lookup(scalar h);
};

#endif
