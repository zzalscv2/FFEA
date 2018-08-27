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
 *      rod_blob_interface.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk
 */
 
#include "rod_blob_interface.h" 

namespace rod{

// interface math

/**
 Given the three nodes in a triangle (or the face of a tetrahedron, if you prefer), get the direction vector of the normal to to that face. 
*/ 
void get_tri_norm(float node0[3], float node1[3], float node2[3], OUT float tri_norm[3]){
    float ax = node1[0] - node0[0];
	float ay = node1[1] - node0[1];
	float az = node1[2] - node0[2];
	float bx = node2[0] - node1[0];
	float by = node2[1] - node1[1];
	float bz = node2[2] - node1[2];
    tri_norm[0] = az*by - ay*bz;
    tri_norm[1] = ax*bz - az*bx;
    tri_norm[2] = ay*bx - ax*by;
}

/**
 For an array of 9 tet_nodes (defined in  object defined in mesh_node.cpp/mesh_node.h) populate a 1-D array with the 3x3 Jacobian of the shape functions. 
*/
void get_jacobian(mesh_node **tet_nodes, OUT float J[0]){
    J[0] = tet_nodes[1]->pos[0] - tet_nodes[0]->pos[0];
    J[1] = tet_nodes[1]->pos[1] - tet_nodes[0]->pos[1];
    J[2] = tet_nodes[1]->pos[2] - tet_nodes[0]->pos[2];
    
    J[3] = tet_nodes[2]->pos[0] - tet_nodes[0]->pos[0];
    J[4] = tet_nodes[2]->pos[1] - tet_nodes[0]->pos[1];
    J[5] = tet_nodes[2]->pos[2] - tet_nodes[0]->pos[2];
    
    J[6] = tet_nodes[3]->pos[0] - tet_nodes[0]->pos[0];
    J[7] = tet_nodes[3]->pos[1] - tet_nodes[0]->pos[1];
    J[8] = tet_nodes[3]->pos[2] - tet_nodes[0]->pos[2];
}

/**
 Find the inverse of a 3x3 matrix (m), and populate the 1-D array m_inv. This is identical to the matrix inversion in mat_vec_fns, only it operates on 1-D arrays of floats instead of 2-D ones of scalars.
*/ 
void float_3x3_invert(float m[9], OUT float m_inv[9]) {
    scalar det;

    // Construct the inverse matrix
    m_inv[0] = m[8] * m[4] - m[7] * m[5];
    m_inv[1] = m[7] * m[2] - m[8] * m[1];
    m_inv[2] = m[5] * m[1] - m[4] * m[2];
    m_inv[3] = m[6] * m[5] - m[8] * m[3];
    m_inv[4] = m[8] * m[0] - m[6] * m[2];
    m_inv[5] = m[3] * m[2] - m[5] * m[0];
    m_inv[6] = m[7] * m[3] - m[6] * m[4];
    m_inv[7] = m[6] * m[1] - m[7] * m[0];
    m_inv[8] = m[4] * m[0] - m[3] * m[1];

    // calc determinant
    det = m[0] * m_inv[0] + m[3] * m_inv[1] + m[6] * m_inv[2];

    // divide by determinant
    det = 1.0 / det;
    m_inv[0] *= det;
    m_inv[1] *= det;
    m_inv[2] *= det;
    m_inv[3] *= det;
    m_inv[4] *= det;
    m_inv[5] *= det;
    m_inv[6] *= det;
    m_inv[7] *= det;
    m_inv[8] *= det;
}

/**
 Multiply the 3x3 matrices A and B. This one is unrolled. Store the result in result. A, B and result are 1-D arrays of length 9 representing 3x3 matrices. Again, this does the same thing as the mat_vec_fn, but with 1-D arrays of floats.
*/ 
void float_3x3_mult_unrolled(float A[9], float B[9], OUT float result[9]){
    result[0] += A[0] * B[0];
    result[0] += A[3] * B[1];
    result[0] += A[6] * B[2];
    result[1] += A[0] * B[3];
    result[1] += A[3] * B[4];
    result[1] += A[6] * B[5];
    result[2] += A[0] * B[6];
    result[2] += A[3] * B[7];
    result[2] += A[6] * B[8];
    result[3] += A[1] * B[0];
    result[3] += A[4] * B[1];
    result[3] += A[7] * B[2];
    result[4] += A[1] * B[3];
    result[4] += A[4] * B[4];
    result[4] += A[7] * B[5];
    result[5] += A[1] * B[6];
    result[5] += A[4] * B[7];
    result[5] += A[7] * B[8];
    result[6] += A[2] * B[0];
    result[6] += A[5] * B[1];
    result[6] += A[8] * B[2];
    result[7] += A[2] * B[3];
    result[7] += A[5] * B[4];
    result[7] += A[8] * B[5];
    result[8] += A[2] * B[6];
    result[8] += A[5] * B[7];
    result[8] += A[8] * B[8];
}

/**
 Transpose a 3x3 matrix, represented as a 1-D array of floats. Populate the array 'transposed'.
*/
void transpose_3x3(float in[9], OUT float transposed[9]){
    transposed[0] = in[0];
    transposed[1] = in[3];
    transposed[2] = in[6];
    transposed[3] = in[1];
    transposed[4] = in[4];
    transposed[5] = in[7];
    transposed[6] = in[2];
    transposed[7] = in[5];
    transposed[8] = in[8];
}

// note: this is mostly the same as the one in tetra_element_linear, but that one is tied to FFEA's data structures,
// objects and types, whereas this one uses pure functions and floats.
/**
 Get the 3x3 gradient deformation matrix of a tetrahedron. The tetrahedron is given as an array of the mesh_node object defined in mesh_node.cpp. Two arrays are needed, one representing the tetrahedron before it was deformed, the other one representing the tetrahedron after deformation.
 \f[ \mathbf{F_e} = (\mathbf{J'} \mathbf{J}^{-1})^T \f]
 Where \f$ \mathbf{F} \f$ is the gradient deformation matrix, \f$ \mathbf{J'} \f$ is the Jacobian of the new deformed tetrahedron, and \f$ \mathbf{J} \f$ is the Jacobian of the undeformed one.
 Note: this is mostly the same as the one in tetra_element_linear.cpp, but that one is a member function that operates on member variables of tetra_element_linear, whereas this one is a pure function (and uses the same floats and 1-d arrays as the rest of the rod stuff).
*/
void get_gradient_deformation(mesh_node **nodes_before, mesh_node**nodes_after, OUT float transposed_gradient_deformation_3x3[9]){
    float J_before[9];
    float J_after[9];
    get_jacobian(nodes_before, J_before);
    get_jacobian(nodes_after, J_after);
    float before_inverse[9];
    float_3x3_invert(J_before, before_inverse);
    float gradient_deformation_3x3[9] = {0,0,0,0,0,0,0,0,0};
    float_3x3_mult_unrolled(J_after, before_inverse, gradient_deformation_3x3);
    transpose_3x3(gradient_deformation_3x3, transposed_gradient_deformation_3x3);
}

/**
 QR decompose\factorise a 3x3 matrix (given here as a 1-d array of floats) using the Gram-Schmidt process. Like the matrix multiplication above, this is again unrolled. This populates arrays for both Q and R, although only R is used.
 We wish to decompose the 3x3 matrix A ([a_1 ... a_n]) into two 3x3 matrices, Q and R.
 The elements e_n of our Q matrix are given by the following recursive definition:
 
 \f[ u_1 = a_1, ~ e_1 = \frac{u1}{||u_1||} \f]
 \f[ u_{k+1} = a_{k+1} - (a_{k+1} \cdot e_1)e1 - ... - (a_{k+1} \cdot e_k)e_k, ~ e_1 = \frac{u1}{||u_1||}, e_{k+1} = \frac{u_{k+1}}{||u_{k+1}||} \f]
 
 \f[A = [ a_1 | a_2 | \dots | a_n] = [ e_1 | e_2 | \dots | e_n]
    \begin{bmatrix} 
    a_{1}\cdot e_1 & a_2\cdot e_1 & \dots & a_n \cdot e_1 \\
    0 & a_2 \cdot e_2 & \dots & a_n \cdot e_2 \\ 
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & a_n \cdot e_n 
    \end{bmatrix}

  = QR \f]
*/ 
void QR_decompose_gram_schmidt(float matrix_3x3[9], OUT float Q[9], float R[9]){
    float a1[3] = {matrix_3x3[0], matrix_3x3[3], matrix_3x3[6]};
    float a2[3] = {matrix_3x3[1], matrix_3x3[4], matrix_3x3[7]};
    float a3[3] = {matrix_3x3[2], matrix_3x3[5], matrix_3x3[8]};
    float u1[3] = {a1[0], a1[1], a1[2]};
    float e1[3];
    normalize(u1, e1);
    float u2[3];
    float a2_dot_e1 = a2[0]*e1[0] + a2[1]*e1[1] + a2[2]*e1[2];
    vec3d(n){ u2[n] = a2[n] - a2_dot_e1*e1[n];}
    float e2[3];
    normalize(u2, e2);
    float u3[3];
    float a3_dot_e1 = a3[0]*e1[0] + a3[1]*e1[1] + a3[2]*e1[2];
    float a3_dot_e2 = a3[0]*e2[0] + a3[1]*e2[1] + a3[2]*e2[2];
    vec3d(n) u3[n] = a3[n] - a3_dot_e1*e1[n] - a3_dot_e2*e2[n];
    float e3[3];
    normalize(u3, e3);
    
    Q[0] = e1[0]; Q[1] = e1[1]; Q[2] = e1[2]; Q[3] = e2[0]; Q[4] = e2[1]; Q[5] = e2[2]; Q[6] = e3[0]; Q[7] = e3[1]; Q[8] = e3[2];
    R[0] = a1[0]*e1[0] + a1[1]*e1[1] + a1[2]*e1[2];
    R[1] = a2[0]*e1[0] + a2[1]*e1[1] + a2[2]*e1[2];
    R[2] = a3[0]*e1[0] + a3[1]*e1[1] + a3[2]*e1[2];
    R[3] = 0;
    R[4] = a2[0]*e2[0] + a2[1]*e2[1] + a2[2]*e2[2];
    R[5] = a3[0]*e2[0] + a3[1]*e2[1] + a3[2]*e2[2];
    R[6] = 0;
    R[7] = 0;
    R[8] = a3[0]*e3[0] + a3[1]*e3[1] + a3[2]*e3[2];
    
    //float Q[9] = { e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2] }
    //float R[9] = { 
    //    a1[0](e1[0] + a1[1](e1[1] + a1[2](e1[2], a2[0](e1[0] + a2[1](e1[1] + a2[2](e1[2], a3[0](e1[0] + a3[1](e1[1] + a3[2](e1[2],
    //    0,  a2[0](e2[0] + a2[1](e2[1] + a2[2](e2[2], a3[0](e2[0] + a3[1](e2[1] + a3[2](e2[2]
    //    0, 0, a3[0](e3[0] + a3[1](e3[1] + a3[2](e3[2]
}

// X1Y2Z3 :)
/**
 Construct an euler rotation matrix of the X1Y2Z3 type. Parameters, a, b and c, the angles of rotation. Populates the array 'rotmat'.
*/
void construct_euler_rotation_matrix(float a, float b, float g, float rotmat[9]){
    float s1 = std::sin(a);
    float s2 = std::sin(b);
    float s3 = std::sin(g);
    float c1 = std::cos(a);
    float c2 = std::cos(b);
    float c3 = std::cos(g);
    rotmat[0] = c2*c3;
    rotmat[1] = -c2*s3;
    rotmat[2] = s2;
    rotmat[3] = c1*s3 + c3*s1*s2;
    rotmat[4] = c1*c3 - s1*s2*s3;
    rotmat[5] = -c2*s1;
    rotmat[6] = s1*s3 - c1*c3*s2;
    rotmat[7] = c3*s1 + c1*s2*s3;
    rotmat[8] = c1*c2;
}

/**
 Rotate a tetrahedron (given by an array of mesh_node objects) by a 3x3 rotation matrix (given as a 1-d array of floats) about its centroid. Populates a new mesh_node object. Note that this function is not used in the main simulation loop, it is only used in the rod_blob_interface unit test.
*/
void rotate_tet(float rotmat[9], mesh_node **nodes, OUT mesh_node **rotated_nodes){
    
    print_array("connected tetrahedron node 0", nodes[0]->pos.data, 3);
    print_array("connected tetrahedron node 1", nodes[1]->pos.data, 3);
    print_array("connected tetrahedron node 2", nodes[2]->pos.data, 3);
    print_array("connected tetrahedron node 3", nodes[3]->pos.data, 3);
    
    // get centroid
    float tet_centroid[3] = {0,0,0};
    for (int i=0; i<4; i++){
        vec3d(n) {tet_centroid[n] += nodes[i]->pos.data[n];}
        
    }
    
    vec3d(n){ tet_centroid[n] /= 4.0; }
    print_array("tet_centroid", tet_centroid, 3);
    
    // move to centre
    for (int i=0; i<4; i++){
        vec3d(n){ nodes[i]->pos[n] -= tet_centroid[n]; }
    }
    
    float posdata_tofloat[3];
    float posdata_rotated[3];

    // rotate
    for (int i=0; i<4; i++){
        print_array("posdata_toflnfr", nodes[i]->pos.data, 3);
        vec3d(n){posdata_tofloat[n] = (float)nodes[i]->pos[n];}
        apply_rotation_matrix(posdata_tofloat, rotmat, posdata_rotated);
        vec3d(n){rotated_nodes[i]->pos.data[n] = posdata_rotated[n];}
        print_array("rotated_nodes[i]->pos.data[n]", rotated_nodes[i]->pos.data, 3);
    }
    
    // move back
    for (int i=0; i<4; i++){
        vec3d(n){ nodes[i]->pos[n] += tet_centroid[n]; }
        vec3d(n){ rotated_nodes[i]->pos[n] += tet_centroid[n]; }
    }

}

// interface structure

/**
 Rod_blob_interface constructor. Requires pointers to fully-initialized rod and blob objects. Parameters:
 - set_ends_at_rod: if true, the connection goes from blob to rod. If false, it goes from rod to blob. This affects the index of the rod element which is connected.
 - set_to_index: index of the source element (via blob->get_element).
 - set_from_index: index of the destination element (via blob->get_element).
 - set_face_index: index of the face of the blob to connected the rod to. Even though we already have the element index (either to 
*/
Rod_blob_interface::Rod_blob_interface(Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int set_face_index)
    {
        this->connected_rod = set_connected_rod;
        this->connected_blob = set_connected_blob;
        this->ends_at_rod = set_ends_at_rod;
        this->to_index = set_to_index;
        this->from_index = set_from_index;
        this->face_index = set_face_index;
        
        if (ends_at_rod){
            this->connected_tet = this->connected_blob->get_element(this->from_index);
        }
        else{
            this->connected_tet = this->connected_blob->get_element(this->to_index);
        }
        
        this->connected_face = this->connected_blob->absolutely_get_face(this->face_index*4); // second order faces are numbered oddly
        
        if (this->connected_face->e != this->connected_tet){
            std::cout << "Face index = " << this->face_index << "\n";
            std::cout << "To index = " << this->to_index << "\n";
            std::cout << "From index = " << this->from_index << "\n";
            rod_abort("The face at the rod_blob_interface does not belong to the tetrahedron of the rod_blob_interface.");
        }
        
        for(int i =0; i<4; i++){
            this->deformed_tet_nodes[i] = new mesh_node();
        }
                
    }

/**
 Update the internal state of the rod-blob interface. The update_edge_vecs parameter will reconstruct the edges, which are need for certain mathematical operations, including updating the position of the attachment node. update_tet updates the cached tetrahedron object stored by the interface.
 
 If you've deformed the tetrahedron somehow, you should leave this alone until you've computed the gradient deformation matrix, as you need 'before' and 'after' tetrahedra to make that work. If you're doing some numerical integraion, you should use the rod_blob_interface tetrahedron as the one you deform, because you can just update_tet as soon as you're done with it.
*/
void Rod_blob_interface::update_internal_state(bool update_edge_vecs, bool update_tet){
    if (update_edge_vecs){this->set_edge_vecs();}
    if (update_tet){this->set_tet(this->connected_tet);}
    std::cout << "Updated internal state. \n";
    print_array("Edge vecs node 0", this->edge_vecs[0], 3);
    print_array("Edge vecs node 1", this->edge_vecs[1], 3);
    print_array("Edge vecs node 2", this->edge_vecs[2], 3);
    print_array("deformed tetrahedron node 0", this->deformed_tet_nodes[0]->pos.data, 3);
    print_array("deformed tetrahedron node 1", this->deformed_tet_nodes[1]->pos.data, 3);
    print_array("deformed tetrahedron node 2", this->deformed_tet_nodes[2]->pos.data, 3);
    print_array("deformed tetrahedron node 3", this->deformed_tet_nodes[3]->pos.data, 3);
}    

/**
 Recompute te edge vectors for a tetrahedron. These are not automatically updated when the nodes are moved. If you do something to the attachment tetrahedron, consider running Rod_blob_interface->update_internal_state().
*/
void Rod_blob_interface::set_edge_vecs(){
    vec3d(v){this->tet_origin[v] = this->connected_tet->n[0]->pos.data[v];}
    vec3d(v){this->edge_vecs[0][v] = this->connected_tet->n[1]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[1][v] = this->connected_tet->n[2]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[2][v] = this->connected_tet->n[3]->pos[v] - this->tet_origin[v];}
}

/**
 Get the position and direction vectors for the current attachment node. You should run this every frame.
*/
void Rod_blob_interface::get_attachment_node(OUT float attachment_node[3], float attachment_node_pos[3]){
    //grab face nodes
    float face_node_1[3];
    float face_node_2[3];
    float face_node_3[3];
    for (int i=0; i<3; i++){
        face_node_1[i] = this->connected_face->n[0]->pos[i];
        face_node_2[i] = this->connected_face->n[1]->pos[i];
        face_node_3[i] = this->connected_face->n[2]->pos[i];
    }
    rod::print_array("face node 1", face_node_1, 3);
    rod::print_array("face node 2", face_node_2, 3);
    rod::print_array("face node 3", face_node_3, 3);
    
    //get normal
    get_tri_norm(face_node_1, face_node_2, face_node_3, attachment_node);
    normalize(attachment_node, attachment_node);
    
    //get position of node in face
    vec3d(n){
        attachment_node_pos[n] = this->tet_origin[n] + (this->edge_vecs[0][n]*this->node_weighting[0] + this->edge_vecs[1][n]*this->node_weighting[1] + this->edge_vecs[2][n]*this->node_weighting[2]);
    }
    
    //make sure normal is facing the right way
    float way_to_centroid[3];
    vec3d(n){way_to_centroid[n] = this->connected_tet->centroid[n] - attachment_node_pos[n];}
    normalize(way_to_centroid, way_to_centroid);
    float dotprod = way_to_centroid[0]*attachment_node[0] + way_to_centroid[1]*attachment_node[1] + way_to_centroid[2]*attachment_node[2];
    if (dotprod < 0){
        vec3d(n){attachment_node[n]*=-1;}
    }
    
}

// need to call this AFTER getting the attachment node
/**
 Once the attachment node has been obtained, this sets the initial direction of the attachment material axis. Note that the initial direction is arbitray - by default it is set such that there is no energy between the attachment axis and the first rod material axis. Only run this function for initialisation! When the simulation is running, you need to update the orientation of this thing using this->reorientate_connection.
*/
void Rod_blob_interface::get_attachment_material_axis(float attachment_node[3], OUT float attachment_material_axis[3]){
    float nearest_material_axis[3];
    float nearest_element[3];
    
    if (this->ends_at_rod){
        vec3d(n){ nearest_material_axis[n] = this->connected_rod->equil_m[this->connected_rod->length - 3 + n]; }
        vec3d(n){ nearest_element[n] = (this->connected_rod->equil_r[this->connected_rod->length - 3 + n]) - (this->connected_rod->equil_r[this->connected_rod->length - 6 + n]); }
    }
    
    else{
        vec3d(n){ nearest_material_axis[n] = this->connected_rod->equil_m[n]; }
        vec3d(n){ nearest_element[n] = this->connected_rod->equil_r[3 + n] - this->connected_rod->equil_r[n]; }
    }
    
    normalize(nearest_material_axis, nearest_material_axis);
    normalize(nearest_element, nearest_element);
    parallel_transport(nearest_material_axis, attachment_material_axis, nearest_element, attachment_node);
    
}

/**
 Overwrite the internal tetrahedron nodes with the nodes from a tetra_element_linear object. 
*/
void Rod_blob_interface::set_tet(tetra_element_linear *tet){
    vec3d(n){this->deformed_tet_nodes[0]->pos[n] = tet->n[0]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[1]->pos[n] = tet->n[1]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[2]->pos[n] = tet->n[2]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[3]->pos[n] = tet->n[3]->pos[n];}
}

/**
 For a given attachment element and material axis, get a new material axis and attachment element based on the rotation matrix found from the QR decomposition of the gradient deformation matrix. Run this every timestep!
*/
void Rod_blob_interface::reorientate_connection(float attachment_element[3], float attachment_material_axis[3], OUT float new_attachment_element[3], float new_attachment_material_axis[3]){
    float gradient_deformation[9] = {0,0,0,0,0,0,0,0,0};
    float R[9];
    float Q[9];
    get_gradient_deformation(this->connected_tet->n, this->deformed_tet_nodes, gradient_deformation);
    QR_decompose_gram_schmidt(gradient_deformation, Q, R);
    apply_rotation_matrix(attachment_element, Q, new_attachment_element); // yes, Q really is the rotation matrix.
    apply_rotation_matrix(attachment_material_axis, Q, new_attachment_material_axis);
}


}
