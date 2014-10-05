#ifndef MATRIXFIXEDSPARSITYPATTERN_H_INCLUDED
#define MATRIXFIXEDSPARSITYPATTERN_H_INCLUDED

#include <vector>

using namespace std;

typedef struct
{
	int x, n;
} sparse_count;

class MatrixFixedSparsityPattern
{
	public:
		int init(tetra_element_linear *elem, int num_elements)
		{
			vector< vector<sparse_count> > all_entries;

			int n, ni, nj;
			for(n = 0; n < num_elements; n++) {
				// add mass matrix for this element
				for(i = 0; i < 4; i++) {
					for(j = 0; j < 4; j++) {
						ni = elem[n].n[i]->index;
						nj = elem[n].n[j]->index;
					}
				}
			}

		}

};

#endif
