#include "mat_vec_types.h"
#include "Face.h"
#include "Blob.h"
#include "mpi.h"
#include <stdlib.h>

struct FaceMpi {
	// each face has a unique index and daddy_blob_index pairs
	int index;
	// int global_index;
	vector3 normal; //
	int daddy_blob_index;
	
	bool kinetically_active;
	vector3 centroid; //

	// arrays with length of num_blobs 
	bool vdw_bb_interaction_flag; //false
	scalar vdw_bb_energy; //0
	vector3 force[3]; 
	vector3 vdw_bb_force; //0
	int n_index[4];
	vector3 n_pos[4];  //
}; 


struct LinkedListNodeMpi {
	FaceMpi *face;
	LinkedListNodeMpi *next;
	int index;
	int x,y,z;
};

class DerivedDatatype {
public:
	 
	// allocate 2d continous array
	scalar **alloc_2d_scalar(int rows, int cols);
	vector3 **alloc_2d_vector3(int rows, int cols);
	FaceMpi **alloc_2d_facempi(int rows, int cols);

	//alternative way of f1 to f2	
	void assign_face(FaceMpi f1, FaceMpi f2);
	
	//create a face from facempi
	void facempi_to_face(FaceMpi *fm, Face *f, int num_blobs);
	
	//create facempi from face
	void face_to_facempi(Face *f, FaceMpi *fm, int num_blobs);
	
	MPI_Datatype derived_scalar_mpi();

	MPI_Datatype derived_vector3_mpi();

	MPI_Datatype derived_face_mpi();
	
};


