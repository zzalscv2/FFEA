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

#ifndef TETRAHEDRAOVERLAP_H_INCLUDED
#define TETRAHEDRAOVERLAP_H_INCLUDED

#ifdef USE_DOUBLE_LESS
typedef float scalar;
#else
typedef double scalar;
#endif

#include <stddef.h>

class checkVars {
	
	public:

		checkVars();
		~checkVars();
		
		// Member variables
		typedef scalar point[3];
		point *V1,*V2;			        ///< vertices coordinates
		scalar e_v1[6][3];            ///< vector edge-oriented 
      scalar e_v2[6][3];            ///< vectors edge-oriented
		int masks[4];  ///< for each face of the first tetrahedron stores the halfspace each vertex of the second tetrahedron belongs to

		scalar P_V1[4][3]; ///< differences between the vertices of the second (first) tetrahedron and the vertex 0  of the first(second) tetrahedron
      scalar P_V2[4][3]; ///< differences between the vertices of the second (first) tetrahedron and the vertex 0  of the first(second) tetrahedron

		scalar  Coord_1[4][4]; ///< vertices coordinates in the affine space
      scalar  Coord_2[4][4]; ///< vertices coordinates in the affine space
		scalar n[3];	  ///< variable to store the normals
};

/** Fast Tetrahedron-Tetrahedron Overlap Algorithm, by Fabio Ganovelli, Frederico Ponchio, Claudio Rocchini. ACM 2002. */
bool tet_a_tet(scalar (&V_1)[4][3],  /* [in] pointers on 3D coord of tetrahedron A */
               scalar (&V_2)[4][3] ); /* [in] pointers on 3D coord of tetrahedron B */


#endif 
