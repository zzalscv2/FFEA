#include "ffea_test.h"

int ffea_test::do_ffea_test(std::string filename){
    int result = 1;
    std::ifstream t(filename);
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::cout << "Test: " << buffer.str();
    
    if (buffer.str().find("connection_test") != std::string::npos ){
        result = ffea_test::connection_test();
    }
    return result;
}

int ffea_test::connection_test(){
    std::cout << "Performing connection test...\n";
    World *world;
    world = new World();
    if(world->init("tet_ascii.1.ffea", 0, 0, 1) == FFEA_ERROR) {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed properly." << endl;
    }
    std::cout << "Test world initialised. \n";
    
    rod::Rod_blob_interface *current_interface = world->rod_blob_interface_array[0];
    
    float attachment_node_pos[3];
    float attachment_node[3];
    float attachment_material_axis[3];
    
    std::cout << "Current blob index in ffea_test: " << world->rod_blob_interface_array[0]->connected_blob->blob_index << "\n";
    
    current_interface->update_internal_state(true, true);
    current_interface->get_attachment_node(attachment_node, attachment_node_pos);
    current_interface->get_attachment_material_axis(attachment_node, attachment_material_axis);
    
    float J[9];
    rod::get_jacobian(current_interface->deformed_tet_nodes, J);
    
    for(int i = 0; i<4; i++){rod::print_array("Tetrahedron node", current_interface->connected_tet->n[i]->pos.data, 3);}
    rod::print_array("attachment node position", attachment_node_pos, 3);
    rod::print_array("attachment node", attachment_node, 3);
    rod::print_array("attachment material axis", attachment_material_axis, 3);
    rod::print_array("Undeformed jacobian", J, 9);



    float rotmat[9];
    float rotation_angle = 0.55;
    rod::construct_euler_rotation_matrix(0, 0, rotation_angle, rotmat);
    rod::rotate_tet(rotmat, current_interface->connected_tet->n, current_interface->deformed_tet_nodes);
    
    rod::print_array("connected tetrahedron node 0", current_interface->connected_tet->n[0]->pos.data, 3);
    rod::print_array("connected tetrahedron node 1", current_interface->connected_tet->n[1]->pos.data, 3);
    rod::print_array("connected tetrahedron node 2", current_interface->connected_tet->n[2]->pos.data, 3);
    rod::print_array("connected tetrahedron node 3", current_interface->connected_tet->n[3]->pos.data, 3);
    
    rod::print_array("rotation matrix", rotmat, 9);
    rod::print_array("Rotated tetrahedron node 1", current_interface->deformed_tet_nodes[0]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 2", current_interface->deformed_tet_nodes[1]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 3", current_interface->deformed_tet_nodes[2]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 4", current_interface->deformed_tet_nodes[3]->pos.data, 3);
    
    float new_attachment_node[3];
    float new_attachment_material_axis[3];
    
    current_interface->reorientate_connection(attachment_node, attachment_material_axis, new_attachment_node, new_attachment_material_axis);
    rod::print_array("Rotated node", new_attachment_node, 3);
    rod::print_array("Rotated material axis", new_attachment_material_axis, 3);
    float dotprod = new_attachment_material_axis[0]*attachment_material_axis[0] + new_attachment_material_axis[1]*attachment_material_axis[1] + new_attachment_material_axis[2]*attachment_material_axis[2];
    float material_axis_rotation_angle = std::acos(dotprod);
    std::cout << "Angle between = " << material_axis_rotation_angle << "\n";

    if ( (material_axis_rotation_angle > rotation_angle -0.01) and (material_axis_rotation_angle < rotation_angle + 0.01)){
        return 0;
    }
    std::cout << "Rod-blob coupling test failed.\n";
    return 1;
}
