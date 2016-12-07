// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

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

    return 0;
}

