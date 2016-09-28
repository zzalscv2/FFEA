#include "DerivedDatatype.h"

scalar** DerivedDatatype::alloc_2d_scalar(int rows, int cols) {
	scalar *data = (scalar *)malloc(rows*cols*sizeof(scalar));
    scalar **array = (scalar **)malloc(rows*sizeof(scalar*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
	
	return array;
}

vector3** DerivedDatatype::alloc_2d_vector3(int rows, int cols) {
	vector3 *data = (vector3 *)malloc(rows*cols*sizeof(vector3));
    vector3 **array = (vector3 **)malloc(rows*sizeof(vector3*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
	
	return array;
}

FaceMpi** DerivedDatatype::alloc_2d_facempi(int rows, int cols) {
	FaceMpi *data = (FaceMpi *)malloc(rows*cols*sizeof(FaceMpi));
    FaceMpi **array = (FaceMpi **)malloc(rows*sizeof(FaceMpi*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
	
	return array;
}

//alternative way of f1 to f2
void DerivedDatatype::assign_face(FaceMpi f1, FaceMpi f2){
	//f1.global_index = f2.global_index;
	f1.index = f2.index;
	f1.normal.x = f2.normal.x;
	f1.normal.y = f2.normal.y;
	f1.normal.z = f2.normal.z;
	f1.daddy_blob_index = f2.daddy_blob_index;
	f1.kinetically_active = f2.kinetically_active;
	f1.centroid.x = f2.centroid.x;
	f1.centroid.y = f2.centroid.y;
	f1.centroid.z = f2.centroid.z;
	for(int i=0;i<3;i++){
		f1.force[i].x = f2.force[i].x;
		f1.force[i].y = f2.force[i].y;
		f1.force[i].z = f2.force[i].z;
	}
	
	f1.vdw_bb_interaction_flag = f2.vdw_bb_interaction_flag;
	f1.vdw_bb_energy = f2.vdw_bb_energy;
	f1.vdw_bb_force.x = f2.vdw_bb_force.x;
	f1.vdw_bb_force.y = f2.vdw_bb_force.y;
	f1.vdw_bb_force.z = f2.vdw_bb_force.z;
	for(int i=0;i<4;i++){
		f1.n_index[i] = f2.n_index[i];
		f1.n_pos[i].x = f2.n_pos[i].x;
		f1.n_pos[i].y = f2.n_pos[i].y;
		f1.n_pos[i].z = f2.n_pos[i].z;
	}
	
}

void DerivedDatatype::facempi_to_face(FaceMpi *fm, Face *f, int num_blobs){
	
		f->index = fm[0].index;
		f->normal.x = fm[0].normal.x;
		f->normal.y = fm[0].normal.y;
		f->normal.z = fm[0].normal.z;
		f->daddy_blob = new Blob();
		f->daddy_blob->blob_index = fm[0].daddy_blob_index;
		f->kinetically_active = fm[0].kinetically_active;
		f->centroid.x = fm[0].centroid.x;
		f->centroid.y = fm[0].centroid.y;
		f->centroid.z = fm[0].centroid.z;
		for(int k=0;k<3;k++){
			f->force[k].x = fm[0].force[k].x;
			f->force[k].y = fm[0].force[k].y;
			f->force[k].z = fm[0].force[k].z; 
		}
		for(int i=0;i<num_blobs;i++){
		f->vdw_bb_interaction_flag[i] = fm[i].vdw_bb_interaction_flag;
		f->vdw_bb_energy[i] = fm[i].vdw_bb_energy;
		f->vdw_bb_force[i].x = fm[i].vdw_bb_force.x;
		f->vdw_bb_force[i].y = fm[i].vdw_bb_force.y;
		f->vdw_bb_force[i].z = fm[i].vdw_bb_force.z;
		for(int k=0;k<4;k++){
			f->n[k] = new mesh_node();
			f->n[k]->index = fm[i].n_index[k];
			f->n[k]->pos.x = fm[i].n_pos[k].x;
			f->n[k]->pos.y = fm[i].n_pos[k].y;
			f->n[k]->pos.z = fm[i].n_pos[k].z;		
		}
	}
}

void DerivedDatatype::face_to_facempi(Face *f, FaceMpi *fm, int num_blobs){
	
		fm[0].index = f->index;
		fm[0].normal.x = f->normal.x;
		fm[0].normal.y = f->normal.y;
		fm[0].normal.z = f->normal.z;
		fm[0].daddy_blob_index = f->daddy_blob->blob_index;
		fm[0].kinetically_active = f->kinetically_active;
		fm[0].centroid.x = f->centroid.x;
		fm[0].centroid.y = f->centroid.y;
		fm[0].centroid.z = f->centroid.z;
		for(int k=0;k<3;k++){
			fm[0].force[k].x = f->force[k].x;
			fm[0].force[k].y = f->force[k].y;
			fm[0].force[k].z = f->force[k].z;
		}
		for(int i=0;i<num_blobs;i++){
		fm[i].vdw_bb_interaction_flag = f->vdw_bb_interaction_flag[i];
		fm[i].vdw_bb_energy = f->vdw_bb_energy[i];
		fm[i].vdw_bb_force.x = f->vdw_bb_force[i].x;
		fm[i].vdw_bb_force.y = f->vdw_bb_force[i].y;
		fm[i].vdw_bb_force.z = f->vdw_bb_force[i].z;
		for(int k=0;k<4;k++){
			fm[i].n_index[k] = f->n[k]->index;
			fm[i].n_pos[k].x = f->n[k]->pos.x;
			fm[i].n_pos[k].y = f->n[k]->pos.y;
			fm[i].n_pos[k].z = f->n[k]->pos.z;		
		}
	}
}

MPI_Datatype DerivedDatatype::derived_scalar_mpi(){
	MPI_Datatype scalar_mpi_type;
	#ifdef USE_DOUBLE
		MPI_Type_contiguous(1, MPI_DOUBLE, &scalar_mpi_type);
	#else
		MPI_Type_contiguous(1, MPI_FLOAT, &scalar_mpi_type);
	#endif
	MPI_Type_commit(&scalar_mpi_type);
	
	return scalar_mpi_type;
}

MPI_Datatype DerivedDatatype::derived_vector3_mpi(){
	//create a mpi derived datatype of vector3
	//one specified field set the upper bound of each element
	MPI_Datatype scalar_mpi_type = derived_scalar_mpi();
	vector3 v;
	int count = 4;
	int	blocklens[4] = {1,1,1,1};
	MPI_Aint displs[4] = {0, (char*)&v.y-(char*)&v,(char*)&v.z-(char*)&v, sizeof(v)};
	MPI_Datatype old_types[4] = {scalar_mpi_type, scalar_mpi_type, scalar_mpi_type, MPI_UB};
	MPI_Datatype vector3_mpi_type;
	MPI_Type_create_struct(count, blocklens, displs, old_types, &vector3_mpi_type);
	MPI_Type_commit(&vector3_mpi_type);
	
	return vector3_mpi_type;
}

MPI_Datatype DerivedDatatype::derived_face_mpi(){
	//create face datatype
	/*face_mpi.vdw_bb_interaction_flag = new bool[params_num_blobs*sizeof(bool)];
	face_mpi.vdw_bb_energy = new scalar[params_num_blobs*sizeof(scalar)];
	face_mpi.vdw_bb_force = new vector3[params_num_blobs*sizeof(vector3)];*/
	FaceMpi f;
	MPI_Datatype scalar_mpi_type = derived_scalar_mpi();
	
	vector3 v;
	MPI_Datatype vector3_mpi_type = derived_vector3_mpi();
	
	int count1 = 12;
	int	blocklens1[12] = {1,1,1,1,1,  1,1,3,1,  4,4,1};
	MPI_Aint displs1[12] = {
		0, 
		(char*)&f.normal-(char*)&f,
		(char*)&f.daddy_blob_index-(char*)&f,
		(char*)&f.kinetically_active-(char*)&f,
		(char*)&f.centroid-(char*)&f,
		(char*)&f.vdw_bb_interaction_flag-(char*)&f,
		(char*)&f.vdw_bb_energy-(char*)&f,
		(char*)&f.force[0]-(char*)&f,
		(char*)&f.vdw_bb_force-(char*)&f,
		(char*)&f.n_index[0]-(char*)&f,
		(char*)&f.n_pos[0]-(char*)&f,
		sizeof(f)
	};
	MPI_Datatype old_types1[12] = {MPI_INT, vector3_mpi_type, MPI_INT, 
									MPI::BOOL, vector3_mpi_type, MPI::BOOL, 
									scalar_mpi_type, vector3_mpi_type, 
									vector3_mpi_type, MPI_INT, vector3_mpi_type,MPI_UB};
	MPI_Datatype face_mpi_type;
	MPI_Type_create_struct(count1, blocklens1, displs1, old_types1, &face_mpi_type);
	MPI_Type_commit(&face_mpi_type);
	
	return face_mpi_type;
}
