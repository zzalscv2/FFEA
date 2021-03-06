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

#ifndef NODEMAP_H_INCLUDED
#define NODEMAP_H_INCLUDED

#include ""
class NodeMap {

	public:
		
		NodeMap(); 

		~NodeMap(); 
	
		int load_map_from_file(int from_conformation, to_conformation, fname);


	// Variables

	int from_conformation; ///< conformation mapped from index

	int to_conformation; ///< conformation mapping to index
};
#endif
