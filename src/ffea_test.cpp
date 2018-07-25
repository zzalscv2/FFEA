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
    
    std::cout << "Current blob index in ffea_test: " << world->rod_blob_interface_array[0]->connected_blob->blob_index << "\n";
    
    current_interface->get_attachment_node(attachment_node, attachment_node_pos);
        
    return 0;
}
