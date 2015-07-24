#include "mat_vec_types.h"
#include <stdlib.h>

class Spring {
public:

    Spring();

    ~Spring();

    /* *  Variables */

    // Spring constant
    scalar k;

    // Equilibrium length
    scalar l;

    // Blobs connected to
    int *blob_index;

    // Conformations connected to
    int *conformation_index;

    // Nodes connected to
    int *node_index;

    // Check if spring is active
    bool am_i_active;
};
