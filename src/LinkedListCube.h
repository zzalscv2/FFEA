#ifndef LINKEDLISTCUBE_H_INCLUDED
#define LINKEDLISTCUBE_H_INCLUDED

#include "FFEA_return_codes.h"

template <typename T>
struct LinkedListNode
{
	/* Pointer to the object this LinkedListNode represents */
	T *obj;

	/* Pointer to the next LinkedListNode in the list */
	LinkedListNode *next;

	/* Remember this node's index in the pool (for matrix construction purposes) */
	int index;

	/* Remember which cell this node has been placed in to avoid recalculations */
	int x, y, z;
};

template <typename T>
class LinkedListCube
{
	public:

		/* Constructor */
		LinkedListCube()
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

		/* Destructor */
		~LinkedListCube()
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

		/* Builds a LinkedListCube of dimensions N_x x N_y x N_z, and an array of LinkedListNodes of size max_num_nodes_in_pool */
		int alloc(int N_x, int N_y, int N_z, int max_num_nodes_in_pool)
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

		/* Adds the specified T object to the pool (at the index given by add_index) */
		int add_to_pool(T *t)
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

		/* Returns pointer to ith object in the pool */
		LinkedListNode<T> * get_from_pool(int i)
		{
			return &pool[i];
		}

		/*
		 * Completely clears entire grid of all linked lists by setting all pointers to NULL in root,
		 * and all 'next' pointers to NULL in the pool of LinkedListNodes.
		 */
		void clear()
		{
			int i;

			// Clear the grid
			for(i = 0; i < N_x*N_y*N_z; i++)
				root[i] = NULL;

			// Clear the pool
			for(i = 0; i < max_num_nodes_in_pool; i++)
				pool[i].next = NULL;
		}

		/*
		 * Adds the specified pool node to the stack at index (x, y, z) on the grid.
		 * Does this by setting the current top of stack node to be the 'next' in the list
		 * after node i, with node i becoming the new top of stack node on that cell;
		 */
		int add_node_to_stack(int i, int x, int y, int z)
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

		/* Returns whatever node is at the top of the (linked list) stack in grid cell (x, y, z) */
		LinkedListNode<T> * get_top_of_stack(int x, int y, int z)
		{
			pbc(&x, &y, &z);

			if(x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
				printf("Error: Looking for stack in out of bounds cell %d %d %d\n", x, y, z);
				return NULL;
			}

			return root[x * N_y * N_z + y * N_z + z];
		}

		/* Returns how many objects are in the 'pool' */
		int get_pool_size()
		{
			return num_nodes_in_pool;
		}

		void get_dim(int *Nx, int *Ny, int *Nz)
		{
			*Nx = N_x;
			*Ny = N_y;
			*Nz = N_z;
		}

	protected:

		/* The cube has dimensions N_x x N_y x N_z */
		int N_x, N_y, N_z;

		/* The number of linked list nodes available for stacking */
		int max_num_nodes_in_pool;

		/* The number of nodes in use */
		int num_nodes_in_pool;

		/* A cubic grid of pointers */
		LinkedListNode<T> **root;

		/* A "pool" of LinkedListNodes available for stacking on the various stacks in the cube */
		LinkedListNode<T> *pool;

	private:
		/* The current index in the pool array at which nodes are to be added */
		int add_index;

		void pbc(int *x, int *y, int *z)
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
};

#endif
