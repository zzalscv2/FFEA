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

void get_tri_norm(float node0[3], float node1[3], float node2[3], OUT float tri_norm[3]);
void get_jacobian(mesh_node **tet_nodes, OUT float J[9]);
void float_3x3_invert(float m[9], OUT float m_inv[9]);
void float_3x3_mult_transpose(float A[9], float B[9], OUT float result[9]);
void get_gradient_deformation(mesh_node **nodes_before, mesh_node**nodes_after, OUT float gradient_deformation_3x3[9]);
void QR_decompose_gram_schmidt(float matrix_3x3[9], OUT float Q[9], float R[9]);
void construct_euler_rotation_matrix(float a, float b, float g, float rotmat[9]);
void rotate_tet(float rotmat[9], mesh_node **nodes, OUT mesh_node **rotated_nodes);

// objects go here

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
    mesh_node* deformed_tet_nodes[4]; // used for calculating the jacobian
    
    //methods
    Rod_blob_interface (Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int set_face_index);
    void update_internal_state(bool update_edge_vecs, bool update_tet);
    void set_edge_vecs();
    void get_attachment_node(OUT float attachment_node[3], float attachment_node_pos[3]);

    void get_initial_material_axis(); // parallel transport nearest mataxis
    void get_attachment_material_axis(float attachment_node[3], OUT float attachment_material_axis[3]);
    void get_tet_rotation_matrix(float tet_points_before[12], float tet_points_after[12]); // calls to get_gradient_deformation and qr_decompose
    
    void make_tet(); // turn nodes into tet elements
    void set_tet(tetra_element_linear *tet);
    void reorientate_connection(float attachment_element[3], float attachment_material_axis[3], OUT float new_attachment_element[3], float new_attachment_material_axis[3]);

    
    // note: maybe a wrapper function for doubles that converts them to floats? see if it works
    
    // note: use vector12 class for internal nodes or something similar
    
};

bool point_inside_tetrahedron(float point[3], float tet[12]); // work out what direction the attachment node should face (it should point inside the tet, i think?)

} //end namespace


#endif
