#ifndef MATRIXFIXEDSPARSITYPATTERN_H_INCLUDED
#define MATRIXFIXEDSPARSITYPATTERN_H_INCLUDED

#include <vector>
#include "tetra_element_linear.h"

using namespace std;

typedef struct
{
	int x, n;
} sparse_count;

class MatrixFixedSparsityPattern
{
	public:
		int init(tetra_element_linear *elem, int num_elements);

};

#endif
