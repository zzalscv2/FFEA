#ifndef ROD_BLOB_INTERFACE
#define ROD_BLOB_INTERFACE

#include "rod_structure.h"

// this bit
#include "SimulationParams.h"
//#include "tetra_element_linear.h"
//#include "Face.h"
//#include "Blob.h"

namespace rod {

// functions go here

struct Rod_blob_interface
{
    // member variables
    Rod* connected_rod;
//    Blob* connected_blob; // commented out for the time being
    bool ends_at_rod = true; // if this is false, it goes rod->blob, otherwise it goes blob->rod
    int to_index;
    int from_index;
    float node_weighting[3] = {0.333333333333333, 0.333333333333333, 0.333333333333333};
    float euler_angles[3] = {0, 0, 0};
    
    //methods
    Rod_blob_interface (Rod* set_connected_rod, /*Blob* set_connected_blob */ bool set_ends_at_rod, int set_to_index, int set_from_index);
};

} //end namespace


#endif
