#ifndef LINKEDLISTCUBE_H_INCLUDED
#define LINKEDLISTCUBE_H_INCLUDED

#include "FFEA_return_codes.h"
#include "Face.h"

template <class T>
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

template <class T>
class LinkedListCube
{
	public:

		/* Constructor */
		LinkedListCube();

		/* Destructor */
		~LinkedListCube();

		/* Builds a LinkedListCube of dimensions N_x x N_y x N_z, and an array of LinkedListNodes of size max_num_nodes_in_pool */
		int alloc(int N_x, int N_y, int N_z, int max_num_nodes_in_pool);

		/* Adds the specified T object to the pool (at the index given by add_index) */
		int add_to_pool(T *t);

		/* Returns pointer to ith object in the pool */
		LinkedListNode<T> * get_from_pool(int i);

		/*
		 * Completely clears entire grid of all linked lists by setting all pointers to NULL in root,
		 * and all 'next' pointers to NULL in the pool of LinkedListNodes.
		 */
		void clear();

		/*
		 * Adds the specified pool node to the stack at index (x, y, z) on the grid.
		 * Does this by setting the current top of stack node to be the 'next' in the list
		 * after node i, with node i becoming the new top of stack node on that cell;
		 */
		int add_node_to_stack(int i, int x, int y, int z);

		/* Returns whatever node is at the top of the (linked list) stack in grid cell (x, y, z) */
		LinkedListNode<T> * get_top_of_stack(int x, int y, int z);

		/* Returns how many objects are in the 'pool' */
		int get_pool_size();

		void get_dim(int *Nx, int *Ny, int *Nz);

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

		void pbc(int *x, int *y, int *z);
};

#include "../src/LinkedListCube.tpp"

#endif
