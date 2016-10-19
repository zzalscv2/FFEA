#ifndef SPARSEMATRIXUNKNOWNPATTERN_H_INCLUDED
#define SPARSEMATRIXUNKNOWNPATTERN_H_INCLUDED

#include <vector>

#include "mat_vec_types.h"
#include "SparseMatrixTypes.h"
#include "FFEA_return_codes.h"

using namespace std;

class SparseMatrixUnknownPattern {
public:
    SparseMatrixUnknownPattern();

    ~SparseMatrixUnknownPattern();

    int init(int num_rows, int suggested_initial_size_for_row_vectors);

    void add_off_diagonal_element(int row_index, int column_index, scalar val);

    void set_diagonal_element(int row_index, scalar val);

    void calc_inverse_diagonal(scalar *inv_D);

    void zero();

    /** Applies this matrix to the given vector 'in', writing the result to 'result' */
    void apply(scalar *in, scalar *result);

    void print();

private:
    int num_rows;
    vector<sparse_entry> *row;
    scalar *diagonal;
};

#endif
