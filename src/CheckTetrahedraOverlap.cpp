/**********************************************************************************
 * tetratetra.cxx                                                        o o      *
 *                                                                o     o         *
 * Visual Computing Group                                         _  O  _         *
 * IEI Institute, CNUCE Institute, CNR Pisa                        \/)\/          *
 *                                                                /\/|            *
 * Copyright(C) 2002 by Fabio Ganovelli, Federico Ponchio and Claudio Rocchini |  *
 * All rights reserved.                                              \            *
 *                                                                                *
 * Permission  to use, copy, modify, distribute  and sell this  software and      *
 * its documentation for any purpose is hereby granted without fee, provided      *
 * that  the above copyright notice appear  in all copies and that both that      *
 * copyright   notice  and  this  permission  notice  appear  in  supporting      *
 * documentation. the author makes  no representations about the suitability      *
 * of this software for any purpose. It is provided  "as is" without express      *
 * or implied warranty.                                                           *
 *                                                                                *
 **********************************************************************************/

/* 
   Fast Tetrahedron-Tetrahedron Overlap Algorithm. 
   
   Authors:
   - Fabio Ganovelli, Istituto di Elaborazione dell'Informazione, National Research Council, Pisa, Italy.
   - Frederico Ponchio, Istituto di Elaborazione dell'Informazione, National Research Council, Pisa, Italy.
   - Claudio Rocchini, Istituto di Elaborazione dell'Informazione, National Research Council, Pisa, Italy.
   
   Abstract: 
   We present an algorithm to test two tetrahedra for overlap. The
   algorithm is based on a dimension reduction technique that allows
   the application of the Separating Axis Theorem, thus avoiding part
   of the computation needed to perform the Separating Axis Test.
*/

#include "CheckTetrahedraOverlap.h"

// ----------- 3D algebraic operators -------------


#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#define VECT(res,a,b) {				\
    res[0] = a[1]*b[2]-b[1]*a[2];		\
    res[1] = -a[0]*b[2]+b[0]*a[2];		\
    res[2] = a[0]*b[1]-b[0]*a[1];		\
  }

#define SUB(res,a,b) {				\
    res[0] = a[0]-b[0];				\
    res[1] = a[1]-b[1];				\
    res[2] = a[2]-b[2];				\
  }

#define SUB_DOT(a,b,c) (			\
			(a[0]-b[0])*c[0]+	\
			(a[1]-b[1])*c[1]+	\
			(a[2]-b[2])*c[2]	\
						)

// Functions from the checkVars object

checkVars::checkVars() {
	V1 = NULL;
	V2 = NULL;
	for(int i = 0; i < 3; ++i) {
		n[i] = 0.0;
		for(int j = 0; j < 6; ++j) {
			e_v1[j][i] = 0.0;
			e_v2[j][i] = 0.0;
		}
		for(int j = 0; j < 4; ++j) {
			masks[j] = 0;
			P_V1[j][i] = 0.0;
			P_V2[j][i] = 0.0;
		}
	}

	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			Coord_1[i][j] = 0.0;
			Coord_2[i][j] = 0.0;
		}
	}
}

checkVars::~checkVars() {
	V1 = NULL;
	V2 = NULL;
	for(int i = 0; i < 3; ++i) {
		n[i] = 0.0;
		for(int j = 0; j < 6; ++j) {
			e_v1[j][i] = 0.0;
			e_v2[j][i] = 0.0;
		}
		for(int j = 0; j < 4; ++j) {
			masks[j] = 0;
			P_V1[j][i] = 0.0;
			P_V2[j][i] = 0.0;
		}
	}

	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			Coord_1[i][j] = 0.0;
			Coord_2[i][j] = 0.0;
		}
	}
}

// FaceA ----------------------------------------------------

inline bool FaceA_1(  scalar * Coord,  int & maskEdges, checkVars vars)
{
  
  maskEdges = 000;
  
  if (( Coord[0] = DOT(vars.P_V1[0] , vars.n)) > 0.0) maskEdges = 001;				
  if (( Coord[1] = DOT(vars.P_V1[1] , vars.n)) > 0.0) maskEdges |= 002;	
  if (( Coord[2] = DOT(vars.P_V1[2] , vars.n)) > 0.0) maskEdges |= 004;	
  if (( Coord[3] = DOT(vars.P_V1[3] , vars.n)) > 0.0) maskEdges |= 010;	
  
  
  return (maskEdges == 017);	// if true it means that all of the vertices are out the halfspace
  
  // defined by this face
  
}
// it is the same as FaceA_1, only the values vars.V2[0]-v_ref are used only for the fourth face

// hence they do not need to be stored

inline bool FaceA_2(scalar * Coord,int & maskEdges, checkVars vars)
{
  maskEdges = 000;
  scalar * v_ref = vars.V1[1];
  
  if (( Coord[0] = SUB_DOT(vars.V2[0],v_ref, vars.n)) > 0.0) maskEdges = 001;	
  if (( Coord[1] = SUB_DOT(vars.V2[1],v_ref, vars.n)) > 0.0) maskEdges |= 002; 
  if (( Coord[2] = SUB_DOT(vars.V2[2],v_ref, vars.n)) > 0.0) maskEdges |= 004; 
  if (( Coord[3] = SUB_DOT(vars.V2[3],v_ref, vars.n)) > 0.0) maskEdges |= 010; 
  
  return (maskEdges == 017);	
}



// FaceB --------------------------------------------------------------

inline bool FaceB_1(checkVars vars)
{
  
  return  ((DOT(vars.P_V2[0] , vars.n) > 0.0) &&				
	   (DOT(vars.P_V2[1] , vars.n) > 0.0) &&	
	   (DOT(vars.P_V2[2] , vars.n) > 0.0) &&	
	   (DOT(vars.P_V2[3] , vars.n) > 0.0));
}

inline bool FaceB_2(checkVars vars)
{
  scalar * v_ref = vars.V2[1];
  return	(( SUB_DOT(vars.V1[0],v_ref , vars.n) > 0.0)  &&
		 ( SUB_DOT(vars.V1[1],v_ref , vars.n) > 0.0)  &&
		 ( SUB_DOT(vars.V1[2],v_ref , vars.n) > 0.0)  &&
		 ( SUB_DOT(vars.V1[3],v_ref , vars.n) > 0.0) );
}


// EdgeA -------------------------------------------------------

inline bool EdgeA(const int & f0 , const int & f1, checkVars vars)
{
  
  scalar * coord_f0 = &vars.Coord_1[f0][0];
  scalar * coord_f1 = &vars.Coord_1[f1][0];
  
  int  maskf0 = vars.masks[f0];
  int  maskf1 = vars.masks[f1];
  
  if( (maskf0 | maskf1) != 017) // if there is a vertex of b 
    
    return false;	      // included in (-,-) return false
  
  
  maskf0 &= (maskf0 ^ maskf1);  // exclude the vertices in (+,+)
  
  maskf1 &= (maskf0 ^ maskf1);
  
  // edge 0: 0--1 
  
  if(  ((maskf0 & 001) &&		// the vertex 0 of b is in (-,+) 
	
	(maskf1 & 002)) &&		// the vertex 1 of b is in (+,-)
       
       ( ((coord_f0[1] * coord_f1[0]) - 
	  (coord_f0[0] * coord_f1[1]))  > 0.0 ) )
    // the edge of b (0,1) intersect (-,-) (see the paper)
    
    return false;   
  
  if(	 ((maskf0 & 002) && (maskf1 & 001)) && ( ((coord_f0[1] * coord_f1[0]) - (coord_f0[0] * coord_f1[1]))  < 0.0 ) )
    return false;   
  
  // edge 1: 0--2 
  
  if(  ((maskf0 & 001) && (maskf1 & 004)) && ( ((coord_f0[2] * coord_f1[0]) - (coord_f0[0] * coord_f1[2]))  > 0.0) )
    return false;	
  
  if(  ((maskf0 & 004) && (maskf1 & 001)) && ( ((coord_f0[2] * coord_f1[0]) - (coord_f0[0] * coord_f1[2]))  < 0.0) )
    return false;	
  
  // edge 2: 0--3 
  
  if(  ((maskf0 & 001) &&(maskf1 & 010)) &&( ((coord_f0[3] * coord_f1[0]) - (coord_f0[0] * coord_f1[3]))  > 0.0) )
    return false;
  
  if(  ((maskf0 & 010) && (maskf1 & 001)) &&( ((coord_f0[3] * coord_f1[0]) - (coord_f0[0] * coord_f1[3]))  < 0.0) )
    return false;	
  
  // edge 3: 1--2 
  
  if(  ((maskf0 & 002) && (maskf1 & 004)) 	&& ( ((coord_f0[2] * coord_f1[1]) - (coord_f0[1] * coord_f1[2]))  > 0.0) )
    return false;
  
  if(  ((maskf0 & 004) && (maskf1 & 002)) 	&& ( ((coord_f0[2] * coord_f1[1]) - (coord_f0[1] * coord_f1[2]))  < 0.0) )
    return false;
  
  
  // edge 4: 1--3 
  
  if(  ((maskf0 & 002) && (maskf1 & 010))  && ( ((coord_f0[3] * coord_f1[1]) - (coord_f0[1] * coord_f1[3]))  > 0.0) )
    return false;
  
  if(  ((maskf0 & 010) && (maskf1 & 002)) 	&& ( ((coord_f0[3] * coord_f1[1]) - (coord_f0[1] * coord_f1[3]))  < 0.0) )
    return false;	
  
  // edge 5: 2--3 
  
  if(  ((maskf0 & 004) && (maskf1 & 010))   && ( ((coord_f0[3] * coord_f1[2]) - (coord_f0[2] * coord_f1[3])) > 0.0) )
    return false;
  
  if(  ((maskf0 & 010) && (maskf1 & 004))   && ( ((coord_f0[3] * coord_f1[2]) - (coord_f0[2] * coord_f1[3])) < 0.0) )
    return false;	
  
  return true;	// there exists a separting plane supported by the edge shared by f0 and f1
  
}

// main function

// bool tet_a_tet(scalar &(V_1)[4][3], scalar &(V_2)[4][3] )
bool tet_a_tet(scalar (&V_1)[4][3],  /* [in] pointers on 3D coord of tetrahedron A */
               scalar (&V_2)[4][3] ) /* [in] pointers on 3D coord of tetrahedron B */
{
  
  // First, we must define the variable object for this call (to ensure thread safety)
  //checkVars vars = new checkVars();
  checkVars vars;
  vars.V1 = V_1;
  vars.V2 = V_2;
  
  SUB(vars.P_V1[0] ,vars.V2[0],vars.V1[0]);	
  SUB(vars.P_V1[1] ,vars.V2[1],vars.V1[0]);	
  SUB(vars.P_V1[2] ,vars.V2[2],vars.V1[0]);	
  SUB(vars.P_V1[3] ,vars.V2[3],vars.V1[0]);	
  
  
  SUB(vars.e_v1[0] , vars.V1[1] , vars.V1[0]);	
  SUB(vars.e_v1[1] , vars.V1[2] , vars.V1[0]);
  

  VECT(vars.n , vars.e_v1[0] ,vars.e_v1[1]);		// find the normal to  face 0
  
  
  if(FaceA_1(&vars.Coord_1[0][0],vars.masks[0], vars))	return false; // if FaceA_1 returns true, it means that a separation plane has been found, and thus returns false, so both tetrahedra don't intersect.
  
  
  SUB(vars.e_v1[2],vars.V1[3],vars.V1[0]);
  VECT(vars.n ,vars.e_v1[2] ,  vars.e_v1[0]);
  
  if(FaceA_1(&vars.Coord_1[1][0], vars.masks[1], vars)) 	return false;		
  
  if(EdgeA(0,1, vars)) return false;	
  
  
  VECT(vars.n,  vars.e_v1[1] , vars.e_v1[2]); 
  
  if(FaceA_1(&vars.Coord_1[2][0], vars.masks[2], vars)) 	return false;	
  
  if(EdgeA(0,2, vars)) return false;	
  if(EdgeA(1,2, vars)) return false;  	
  
  SUB(vars.e_v1[4], vars.V1[3],vars.V1[1]);
  SUB(vars.e_v1[3], vars.V1[2],vars.V1[1]);
  
  VECT(vars.n ,vars.e_v1[4] , vars.e_v1[3]);
  
  if(FaceA_2(&vars.Coord_1[3][0],vars.masks[3], vars))  return false;	
  
  if(EdgeA(0,3, vars)) return false;	
  if(EdgeA(1,3, vars)) return false; 	
  if(EdgeA(2,3, vars)) return false; 	
  
  if( (vars.masks[0] | vars.masks[1] | vars.masks[2] | vars.masks[3] )!=017) return true; 
  
  
  // from now on, if there is a separating plane it is parallel to a face of b
  
  SUB(vars.P_V2[0] , vars.V1[0],vars.V2[0]);
  SUB(vars.P_V2[1] , vars.V1[1],vars.V2[0]);	
  SUB(vars.P_V2[2] , vars.V1[2],vars.V2[0]);	
  SUB(vars.P_V2[3] , vars.V1[3],vars.V2[0]);	
  
  
  SUB(vars.e_v2[0] , vars.V2[1], vars.V2[0]);
  SUB(vars.e_v2[1] , vars.V2[2], vars.V2[0]);
  
  VECT(vars.n, vars.e_v2[0] , vars.e_v2[1] );
  if(FaceB_1(vars)) return false;	
  
  SUB(vars.e_v2[2], vars.V2[3], vars.V2[0]);
  
  VECT(vars.n,  vars.e_v2[2] ,  vars.e_v2[0]);
  
  if(FaceB_1(vars)) return false;	
  
  VECT(vars.n,  vars.e_v2[1] ,vars.e_v2[2]);
  
  if(FaceB_1(vars)) return false;
  
  SUB(vars.e_v2[4] , vars.V2[3] , vars.V2[1]);
  SUB(vars.e_v2[3] , vars.V2[2] , vars.V2[1]);
  
  VECT(vars.n , vars.e_v2[4] , vars.e_v2[3]);
  
  if(FaceB_2(vars)) return false;

  return true;	
}

/* *********************************************************** */


int intersectTetrahedronGanovelli(scalar * coordsA[4], 
				  scalar * coordsB[4]) 
{
  scalar V_1[4][3];
  scalar V_2[4][3];
  for(int n = 0 ; n < 4 ; n++)
    for(int d = 0 ; d < 3 ; d++) {
      V_1[n][d] = coordsA[n][d];
      V_2[n][d] = coordsB[n][d];
    }
  
  if(tet_a_tet(V_1,V_2)) 
    return 1;
  
  return 0;
}

/* *********************************************************** */
