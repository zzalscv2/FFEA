#ifndef TETRAHEDRAOVERLAP_H_INCLUDED
#define TETRAHEDRAOVERLAP_H_INCLUDED

#ifdef USE_DOUBLE
typedef double scalar;
#else
typedef float scalar;
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
