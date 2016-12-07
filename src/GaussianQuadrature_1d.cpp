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

#include "GaussianQuadrature_1d.h"

scalar GaussianQuadrature_1d::integrate_function_1d(scalar a, scalar b, int precision) {
    scalar cp = (b + a) / 2.0, cm = (b - a) / 2.0;
    scalar I = 0;
    for (int i = gq_precision_index[precision]; i < gq_precision_index[precision] + precision; i++) {
        I += gq_line[i].W * f_1d(cm * gq_line[i].zeta + cp);
    }
    return I * cm;
}

scalar GaussianQuadrature_1d::integrate_function_1d_tri(scalar theta_max, scalar L_perp, scalar theta_star, int precision) {
    scalar c = theta_max / 2.0;
    scalar I = 0;
    for (int i = gq_precision_index[precision]; i < gq_precision_index[precision] + precision; i++) {
        I += gq_line[i].W * integrate_function_1d(0, L_perp / (cos((c * (gq_line[i].zeta + 1)) - theta_star)), precision);
    }
    return I * c;
}

