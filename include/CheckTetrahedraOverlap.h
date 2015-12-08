#ifndef TETRAHEDRAOVERLAP_H_INCLUDED
#define TETRAHEDRAOVERLAP_H_INCLUDED

#ifdef USE_DOUBLE
typedef double scalar;
#else
typedef float scalar;
#endif


/** Fast Tetrahedron-Tetrahedron Overlap Algorithm, by Fabio Ganovelli, Frederico Ponchio, Claudio Rocchini. ACM 2002. */
bool tet_a_tet(scalar (&V_1)[4][3],  /* [in] pointers on 3D coord of tetrahedron A */
               scalar (&V_2)[4][3] ); /* [in] pointers on 3D coord of tetrahedron B */


#endif 
