#include "Spring.h"

Spring::Spring() {
    blob_index = new int[2];
    conformation_index = new int[2];
    node_index = new int[2];
    k = 0.0;
    l = 0.0;
    am_i_active = false;
}

Spring::~Spring() {
    delete[] blob_index;
    blob_index = NULL;
    delete[] conformation_index;
    conformation_index = NULL;
    delete[] node_index;
    node_index = NULL;
    k = 0.0;
    l = 0.0;
    am_i_active = false;
}

