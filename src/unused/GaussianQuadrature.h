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

#ifndef GAUSSIANQUADRATURE_H_INCLUDED
#define GAUSSIANQUADRATURE_H_INCLUDED

#include "Face.h"

typedef struct
{
	int index, num_points;
} precision_lookup;

typedef struct
{
	/* Weight */
	scalar W;

	/* Barycentric coords of point on triangle */
	scalar zeta_0, zeta_1, zeta_2;
} barycentric_gq_point;


static barycentric_gq_point gq_triangle[] = {
	/* 1 point - precision 1 */
	{1.0,	1.0/3.0, 1.0/3.0, 1.0/3.0},

	/* 3 point - precision 2 */
	{1.0/3.0,	2.0/3.0, 1.0/6.0, 1.0/6.0},
	{1.0/3.0,	1.0/6.0, 2.0/3.0, 1.0/6.0},
	{1.0/3.0,	1.0/6.0, 1.0/6.0, 2.0/3.0},

	/* 4 point - precision 3 */
	{-0.562500000000000,	1.0/3.0, 1.0/3.0, 1.0/3.0},

	{0.520833333333333,	.6, .2, .2},
	{0.520833333333333,	.2, .6, .2},
	{0.520833333333333,	.2, .2, .6},

	/* 6 point - precision 4 */
	{0.109951743655322,	0.816847572980459, 0.091576213509771, 0.091576213509771},
	{0.109951743655322,	0.091576213509771, 0.816847572980459, 0.091576213509771},
	{0.109951743655322,	0.091576213509771, 0.091576213509771, 0.816847572980459},

	{0.223381589678011,	0.108103018168070, 0.445948490915965, 0.445948490915965},
	{0.223381589678011,	0.445948490915965, 0.108103018168070, 0.445948490915965},
	{0.223381589678011,	0.445948490915965, 0.445948490915965, 0.108103018168070},

	/* 7 point - precision 5 */
	{0.225000000000000,	1.0/3.0, 1.0/3.0, 1.0/3.0},

	{0.125939180544827,	0.797426985353087, 0.101286507323456, 0.101286507323456},
	{0.125939180544827,	0.101286507323456, 0.797426985353087, 0.101286507323456},
	{0.125939180544827,	0.101286507323456, 0.101286507323456, 0.797426985353087},

	{0.132394152788506,	0.059715871789770, 0.470142064105115, 0.470142064105115},
	{0.132394152788506,	0.470142064105115, 0.059715871789770, 0.470142064105115},
	{0.132394152788506,	0.470142064105115, 0.470142064105115, 0.059715871789770},

	/* 12 point - precision 6 */
	{0.050844906370207,	0.873821971016996, 0.063089014491502, 0.063089014491502},
	{0.050844906370207,	0.063089014491502, 0.873821971016996, 0.063089014491502},
	{0.050844906370207,	0.063089014491502, 0.063089014491502, 0.873821971016996},

	{0.116786275726379,	0.501426509658179, 0.249286745170910, 0.249286745170910},
	{0.116786275726379,	0.249286745170910, 0.501426509658179, 0.249286745170910},
	{0.116786275726379,	0.249286745170910, 0.249286745170910, 0.501426509658179},

	{0.082851075618374,	0.636502499121399, 0.310352451033785, 0.053145049844816},
	{0.082851075618374,	0.310352451033785, 0.053145049844816, 0.636502499121399},
	{0.082851075618374,	0.053145049844816, 0.636502499121399, 0.310352451033785},
	{0.082851075618374,	0.636502499121399, 0.053145049844816, 0.310352451033785},
	{0.082851075618374,	0.310352451033785, 0.636502499121399, 0.053145049844816},
	{0.082851075618374,	0.053145049844816, 0.310352451033785, 0.636502499121399}
};

static precision_lookup gq_precision[] = {
	{0,0},	// Ignore this line (no precision 0)
	{0, 1},	// Precision 1: Starts at index 0, 1 point
	{1, 3},	// Precision 2: Starts at index 1, 3 points
	{4, 4},	// Precision 3: Starts at index 4, 4 points
	{8, 6},	// Precision 4: Starts at index 8, 6 points
	{14, 7},// Precision 5: Starts at index 14, 7 points
	{21, 12}// Precision 6: Starts at index 21, 12 points
};

class GaussianQuadrature
{
	public:

		virtual ~GaussianQuadrature() {}

		/* Integrates function f(p,q), for fixed point p and face coordinate q, at the given precision */
		scalar integrate_point_to_face(vector3 *p, Face *face, int precision)
		{
			vector3 q;
			scalar result = 0;
			int j = gq_precision[precision].index;
			for(int i = 0; i < gq_precision[precision].num_points; i++) {
				face->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, &q);
				result += gq_triangle[j].W * f(p, &q);
				j++;
			}
			return (result * face->area);
		}

		/* Integrates function f(p, q) between two faces, at the given precision */
		scalar integrate_face_to_face(Face *f1, Face *f2, int precision)
		{
			vector3 q;
			scalar result = 0;
			int j = gq_precision[precision].index;
			for(int i = 0; i < gq_precision[precision].num_points; i++) {
				f1->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, &q);
				result += gq_triangle[j].W * integrate_point_to_face(&q, f2, precision);
				j++;
			}
			return (result * f1->area);
		}

	protected:

//		static barycentric_gq_point gq_triangle[NUM_GAUSS_QUAD_POINTS];
//		static precision_lookup gq_precision[NUM_PRECISIONS + 1];

		/* The function being integrated (to be defined in implementing sub class) */
		virtual scalar f(vector3 *p, vector3 *q) = 0;
};

#endif
