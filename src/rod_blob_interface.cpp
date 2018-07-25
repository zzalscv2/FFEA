#include "rod_blob_interface.h" 

namespace rod{
    
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
        
        this->connected_face = this->connected_blob->absolutely_get_face(this->face_index);
        
        if (this->connected_face->e != this->connected_tet){
            std::cout << "Face index = " << this->face_index << "\n";
            std::cout << "To index = " << this->to_index << "\n";
            std::cout << "From index = " << this->from_index << "\n";
            rod_abort("The face at the rod_blob_interface does not belong to the tetrahedron of the rod_blob_interface.");
        }
        
        vec3d(v){this->tet_origin[v] = this->connected_tet->n[0]->pos[v];}
        
        this->set_edge_vecs();
                
    }
    
void Rod_blob_interface::set_edge_vecs(){
    vec3d(v){this->edge_vecs[0][v] = this->connected_tet->n[0]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[0][v] = this->connected_tet->n[0]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[0][v] = this->connected_tet->n[0]->pos[v] - this->tet_origin[v];}
}

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
    
    //get normal
    get_tri_norm(face_node_1, face_node_2, face_node_3, attachment_node);
    normalize(attachment_node, attachment_node);
    
    //get position of node in face
    
    vec3d(n){attachment_node[n] = this->tet_origin[n] + (this->edge_vecs[0][n]*node_weighting[0] + this->edge_vecs[1][n]*node_weighting[1] + this->edge_vecs[2][n]*node_weighting[2]);}
    
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
    
    parallel_transport(nearest_material_axis, attachment_material_axis, nearest_element, attachment_node);
    
}






}
