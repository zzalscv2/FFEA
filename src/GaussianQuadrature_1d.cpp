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

