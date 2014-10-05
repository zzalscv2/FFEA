#ifndef SPARSITYPATTERN_H_INCLUDED
#define SPARSITYPATTERN_H_INCLUDED

#include <list>
#include <vector>

#include "mat_vec_types.h"
#include "FFEA_return_codes.h"

#include "SparseMatrixTypes.h"
#include "SparseMatrixFixedPattern.h"

using namespace std;

typedef struct
{
	int column_index;
	vector<scalar *> source_list;
} sparse_contribution_location;

class SparsityPattern
{
	public:
		SparsityPattern()
		{
			num_rows = 0;
			row = NULL;
			num_nonzero_elements = 0;
		}

		~SparsityPattern()
		{
			delete[] row;
			num_rows = 0;
			row = NULL;
			num_nonzero_elements = 0;
		}

		int init(int num_rows)
		{
			this->num_rows= num_rows;
			row = new list<sparse_contribution_location*>[num_rows];
			if(row == NULL) {
				printf("Could not allocate memory for 'row' array in SparsityPattern\n");
				return FFEA_ERROR;
			}

			return FFEA_OK;
		}

		/*
		 *
		 */
		void register_contribution(int i, int j, scalar *contrib_memory_loc)
		{
			list<sparse_contribution_location*>::iterator it;
			for(it = row[i].begin(); it != row[i].end(); it++) {

				// If element already has sources, add the source to the list
				if((*it)->column_index == j) {
					(*it)->source_list.push_back(contrib_memory_loc);
					return;
				}

				// If we've passed the point where we'd expect our element to be,
				// create a new contribution list and insert it
				if((*it)->column_index > j) {
					break;
				}
			}

			num_nonzero_elements++;
			sparse_contribution_location *scl = new sparse_contribution_location();
			scl->column_index = j;
			scl->source_list.push_back(contrib_memory_loc);
			row[i].insert(it, scl);

		}

		bool check_for_contribution(int i, int j) 
		{
			list<sparse_contribution_location*>::iterator it;
			for(it = row[i].begin(); it != row[i].end(); it++) {
				// If element already has sources, add the source to the list
				if((*it)->column_index == j) {
					return true;
				}
			}
			return false;
		}

		/* Factory function for making empty fixed sparsity pattern matrices from this sparsity pattern */
		SparseMatrixFixedPattern * create_sparse_matrix()
		{
			SparseMatrixFixedPattern *sm = new SparseMatrixFixedPattern();

			// Generate the array of sparse matrix elements from this sparsity pattern
			sparse_entry *entry = new sparse_entry[num_nonzero_elements];

			// Also generate the key (array of pointers that take us to the start of each row)
			int *key = new int[num_rows + 1];

			// Initialise it to zero (to get rid of compiler warnings)
			for(int i = 0; i < num_rows; i++) {
				key[i] = 0;
			}

			// Generate the array of sources for each element in the sparse matrix
			sparse_entry_sources *source_list = new sparse_entry_sources[num_nonzero_elements];

			int pos = 0;
			for(int i = 0; i < num_rows; i++) {

				// Store the index of the start of each row in the key
				key[i] = pos;

				// Dump the data serially, initialising values to zero and column indices to those
				// given by the sparsity pattern
				for(list<sparse_contribution_location*>::iterator it = row[i].begin(); it != row[i].end(); it++) {
					entry[pos].val = 0;
					entry[pos].column_index = (*it)->column_index;

					source_list[pos].init((*it)->source_list.size());
					for(int j = 0; j < source_list[pos].get_num_sources(); j++) {
						source_list[pos].set_source(j, (*it)->source_list[j]);
					}

					pos++;
				}
			}

			// Last index in the key array should contain the total number of nonzero elements in the matrix
			key[num_rows] = pos;

			// Initialise the sparse matrix with the array of entries and the row access key
			sm->init(num_rows, num_nonzero_elements, entry, key, source_list);

			// Return pointer to the newly allocated and initialised sparse matrix
			return sm;
		}

		void print()
		{
			for(int i = 0; i < num_rows; i++) {
				printf("= ");
				for(list<sparse_contribution_location*>::iterator it = row[i].begin(); it != row[i].end(); it++) {
					printf("[%d %d] ", (*it)->column_index, (int)((*it)->source_list.size()));
				}
				printf("\n");
			}
		}

	private:

		/* Number of rows in matrix */
		int num_rows;

		/* An array of vectors containing the indices of the occupied sites */
		list<sparse_contribution_location*> *row;

		/* Total number of nonzero elements in the sparsity pattern */
		int num_nonzero_elements;
};

#endif
