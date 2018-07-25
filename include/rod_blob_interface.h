#ifndef ROD_BLOB_INTERFACE
#define ROD_BLOB_INTERFACE

#include "rod_structure.h"

// this bit
//#include "SimulationParams.h"
#include "tetra_element_linear.h"
#include "Face.h"
#include "Blob.h"

namespace rod {

// functions go here

struct Rod_blob_interface
{
    // member variables
    Rod* connected_rod;
    Blob* connected_blob;
    Face* connected_face;
    tetra_element_linear* connected_tet;
    bool ends_at_rod = true; // if this is false, it goes rod->blob, otherwise it goes blob->rod
    int to_index;
    int from_index;
    int face_index;
    float node_weighting[3] = {0.333333333333333, 0.333333333333333, 0.333333333333333};
    float euler_angles[3] = {0, 0, 0};
    float tet_origin[3];
    float edge_vecs[3][3];
    
    
    //methods
    Rod_blob_interface (Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int set_face_index);
    void get_attachment_material_axis(float attachment_node[3], OUT float attachment_material_axis[3]);

    void set_edge_vecs();
    void get_attachment_node(OUT float attachment_node[3], float attachment_node_pos[3]);

    void get_initial_material_axis(); // parallel transport nearest mataxis
    
    void get_tet_rotation_matrix(float tet_points_before[12], float tet_points_after[12]); // calls to get_gradient_deformation and qr_decompose
    
    void make_tet(); // turn nodes into tet elements
    
    // note: maybe a wrapper function for doubles that converts them to floats? see if it works
    
    // note: use vector12 class for internal nodes or something similar
    
};

bool point_inside_tetrahedron(float point[3], float tet[12]); // work out what direction the attachment node should face (it should point inside the tet, i think?)

} //end namespace


#endif
