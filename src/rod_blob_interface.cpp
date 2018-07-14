#include "rod_blob_interface.h" 



namespace rod{

Rod_blob_interface::Rod_blob_interface(Rod* set_connected_rod, /*Blob* set_connected_blob, */ bool set_ends_at_rod, int set_to_index, int set_from_index)
    {
        this->connected_rod = set_connected_rod;
        //this->connected_blob = set_connected_blob;
        this->ends_at_rod = set_ends_at_rod;
        this->to_index = set_to_index;
        this->from_index = set_from_index;
    }; 

}


