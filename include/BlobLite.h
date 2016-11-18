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
  /** BlobLite constructor:
   * Initialises all variables and pointers to 0 (or NULL). 
   * Does not perform any memory allocation. */
  BlobLite();  
  /** BlobLite destructor:
   * Deallocates memory held by arrays, */
  ~BlobLite();
 
  /** Opens and reads the given 'ffea node file', extracting all the nodes for this Blob.
   * Records how many of these are surface nodes and how many are interior nodes.  */
  int load_nodes(const char *node_filename, scalar scale); 
  /** Opens and reads the given 'ffea topology file', extracting all the elements for this Blob.
    * Records how many of these are surface elements and how many are interior elements.  */
  int load_topology(const char *topology_filename);
  /** Read the following num_nodes lines within the given trj file 
    *  and set the coordinates of the nodes accordingly */ 
  int read_nodes_from_file(FILE *trj); 
  /** Calculate and return the center of coordinates */
  int center_of_coord(arr3 &cm); 

  int num_nodes; ///< number of nodes
  int num_surface_nodes; ///< number of surface nodes 
  int num_interior_nodes; ///< number of interior nodes 
  int num_elements; ///< number of elements 
  int num_surface_elements; ///< number of surface elements
  int num_interior_elements; ///< number of interior elements 
  int blob_state; ///< state of the blob, ACTIVE is default.

  scalar *coord; ///< pointer with all nodes coordinates xyzxyzxyz
  int *elem; ///< pointer with all the element indices,  n0i0 n0i1 ... n0i9 n1i0 n1i1 ... n1i9 ... nNi9.

private: 
  /** Part of load_topology, it stores the index of a node within the *elem list */
  int store_index_to_elemnode(int node_index, int elem_number, int node_number); 
  /** Returns the coord index for a node nodei in element elemi to be used in coord[ ] */ 
  int icoord_for_elem_node(int elemi, int nodei); 

};

