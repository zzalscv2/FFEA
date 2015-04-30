#ifndef SPARSEMATRIXTYPES_H_INCLUDED
#define SPARSEMATRIXTYPES_H_INCLUDED

#include "mat_vec_types.h"
#include "FFEA_return_codes.h"

typedef struct {
    /* The column index of this entry in the original matrix */
    int column_index;

    /* The value of this entry */
    scalar val;
} sparse_entry;

class sparse_entry_sources {
public:
    sparse_entry_sources();

    ~sparse_entry_sources();

    int init(int num_sources);

    void set_source(int i, scalar *s);

    scalar sum_all_sources();

    int get_num_sources();

private:
    int num_sources;
    scalar **sources;
};

#endif
