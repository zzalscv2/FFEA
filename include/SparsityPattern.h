#ifndef SPARSITYPATTERN_H_INCLUDED
#define SPARSITYPATTERN_H_INCLUDED

#include <list>
#include <vector>

#include "mat_vec_types.h"
#include "FFEA_return_codes.h"

#include "SparseMatrixTypes.h"
#include "SparseMatrixFixedPattern.h"

using namespace std;

typedef struct {
    int column_index;
    vector<scalar *> source_list;
} sparse_contribution_location;

class SparsityPattern {
public:
    SparsityPattern();

    ~SparsityPattern();

    int init(int num_rows);

    /* * */
    void register_contribution(int i, int j, scalar *contrib_memory_loc);

    bool check_for_contribution(int i, int j);

    /* Factory function for making empty fixed sparsity pattern matrices from this sparsity pattern */
    SparseMatrixFixedPattern * create_sparse_matrix();

    void print();

private:

    /* Number of rows in matrix */
    int num_rows;

    /* An array of vectors containing the indices of the occupied sites */
    list<sparse_contribution_location*> *row;

    /* Total number of nonzero elements in the sparsity pattern */
    int num_nonzero_elements;
};

#endif
