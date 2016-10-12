#include <cstring>
#include <iostream>
#include "mat_vec_types.h"
#include "mat_vec_fns_II.h"
#include "FFEA_return_codes.h"
#include "dimensions.h"

using namespace std;

#define NUM_NODES_QUADRATIC_TET 10

class BlobLite {
public:  
  /**
   * BlobLite constructor:
   * Initialises all variables and pointers to 0 (or NULL). 
   * Does not perform any memory allocation. */
  BlobLite();  
  /**
   * BlobLite destructor:
   * Deallocates memory held by arrays, */
  ~BlobLite();
 
  int load_nodes(const char *node_filename, scalar scale); 
  int load_topology(const char *topology_filename);
  int read_nodes_from_file(FILE *trj); 
  int store_index_to_elemnode(int node_index, int elem_number, int node_number); 
  int icoord_for_elem_node(int elemi, int nodei);
  int num_nodes; 
  int num_surface_nodes; 
  int num_interior_nodes; 
  int num_elements; 
  int num_surface_elements; 
  int num_interior_elements; 
  int blob_state; 
  int center_of_coord(arr3 &cm);

  scalar *coord; // xyzxyzxyz
  int *elem; // n0i0 n0i1 ... n0i9 n1i0 n1i1 ... n1i9 ... nNi9.

};

