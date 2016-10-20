#include "mat_vec_types.h"
#include <stdlib.h>

class Spring {
public:

    Spring();

    ~Spring();

    /* *  Variables */

    scalar k; ///< Spring constant

    scalar l; ///< Equilibrium length

    int *blob_index; ///< Blobs connected to

    int *conformation_index; ///< Conformations connected to

    int *node_index; ///< Nodes connected to

    bool am_i_active; ///< Check if spring is active
};
