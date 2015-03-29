#include "LinkedListCube.h"


/*  */
template <class T> 
LinkedListCube<T>::LinkedListCube()
{
	N_x = 0;
	N_y = 0;
	N_z = 0;
	max_num_nodes_in_pool = 0;
	num_nodes_in_pool = 0;
	add_index = 0;
	root = NULL;
	pool = NULL;
}

		/*  */
template <class T> 
LinkedListCube<T>::~LinkedListCube()
{
	delete[] pool;
	delete[] root;

	N_x = 0;
	N_y = 0;
	N_z = 0;
	max_num_nodes_in_pool = 0;
	num_nodes_in_pool = 0;
	add_index = 0;
	root = NULL;
	pool = NULL;
}

/* */
template <class T> 
int LinkedListCube<T>::alloc(int N_x, int N_y, int N_z, int max_num_nodes_in_pool)
{
	this->N_x = N_x;
	this->N_y = N_y;
	this->N_z = N_z;
	this->max_num_nodes_in_pool = max_num_nodes_in_pool;
	num_nodes_in_pool = 0;
	root = new LinkedListNode<T> * [N_x * N_y * N_z];
	pool = new LinkedListNode<T>[max_num_nodes_in_pool];

	if(root == NULL || pool == NULL) {
		FFEA_ERROR_MESSG("Could not allocate memory (for root and pool arrays) in LinkedListCube\n");
	}

	// Make sure all pointers are initialised to NULL
	clear();

	// Set the current addition index to the beginning of the pool array
	add_index = 0;

	return FFEA_OK;
}

/* */
template <class T> 
int LinkedListCube<T>::add_to_pool(T *t)
{
	if(add_index >= max_num_nodes_in_pool) {
		FFEA_ERROR_MESSG("In LinkedListCube, attempt to add more nodes to the pool than space has been allocated for.\n");
	}

	// Give object its own representant LinkedListNode in the pool
	pool[add_index].obj = t;
	pool[add_index].index = add_index;
	add_index++;

	num_nodes_in_pool++;

	return FFEA_OK;
}

/* */
template <class T> 
LinkedListNode<T> * LinkedListCube<T>::get_from_pool(int i)
{
	return &pool[i];
}

/* */
template <class T> 
void LinkedListCube<T>::clear()
{
	int i;

	// Clear the grid
	for(i = 0; i < N_x*N_y*N_z; i++)
		root[i] = NULL;

	// Clear the pool
	for(i = 0; i < max_num_nodes_in_pool; i++)
		pool[i].next = NULL;
}

/* */
template <class T> 
int LinkedListCube<T>::add_node_to_stack(int i, int x, int y, int z)
{
	// Apply PBC to the face centroid

	pbc(&x, &y, &z);

	if(x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
		printf("Face centroid out of bounds of LinkedListCube (coords are [%d, %d, %d])\n", x, y, z);
		return FFEA_ERROR;
	}

	int abs_index = x * N_y * N_z + y * N_z + z;

	pool[i].x = x;
	pool[i].y = y;
	pool[i].z = z;
	pool[i].next = root[abs_index];
	root[abs_index] = &pool[i];

	return FFEA_OK;
}

/* */
template <class T> 
LinkedListNode<T> * LinkedListCube<T>::get_top_of_stack(int x, int y, int z)
{
	pbc(&x, &y, &z);

	if(x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
		printf("Error: Looking for stack in out of bounds cell %d %d %d\n", x, y, z);
		return NULL;
	}

	return root[x * N_y * N_z + y * N_z + z];
}

/* */
template <class T> 
int LinkedListCube<T>::get_pool_size()
{
	return num_nodes_in_pool;
}

template <class T> 
void LinkedListCube<T>::get_dim(int *Nx, int *Ny, int *Nz)
{
	*Nx = N_x;
	*Ny = N_y;
	*Nz = N_z;
}


template <class T> 
void LinkedListCube<T>::pbc(int *x, int *y, int *z)
{
	if(*x < 0) {
		*x += N_x;
	} else if(*x >= N_x) {
		*x -= N_x;
	}
	if(*y < 0) {
		*y += N_y;
	} else if(*y >= N_y) {
		*y -= N_y;
	}
	if(*z < 0) {
		*z += N_z;
	} else if(*z >= N_z) {
		*z -= N_z;
	}
}

