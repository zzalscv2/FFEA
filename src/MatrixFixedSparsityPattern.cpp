#include "MatrixFixedSparsityPattern.h"

int MatrixFixedSparsityPattern::init(tetra_element_linear *elem, int num_elements) {
    vector< vector<sparse_count> > all_entries;

    int n, ni, nj;
    for (n = 0; n < num_elements; n++) {
        // add mass matrix for this element
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                ni = elem[n].n[i]->index;
                nj = elem[n].n[j]->index;
            }
        }
    }

}

