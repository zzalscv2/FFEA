#include "GaussianQuadrature_tri.h"

/* */
scalar GaussianQuadrature_tri::integrate_point_to_face(scalar(*f)(vector3*, vector3*), vector3 *p, Face *face, int precision) {
    vector3 q;
    scalar result = 0;
    int j = gq_precision[precision].index;
    for (int i = 0; i < gq_precision[precision].num_points; i++) {
        face->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, &q);
        result += gq_triangle[j].W * f_3d(p, &q);
        j++;
    }
    return (result * face->area);
}

/* */
scalar GaussianQuadrature_tri::integrate_face_to_face(scalar(*f)(vector3*, vector3*), Face *f1, Face *f2, int precision) {
    vector3 q;
    scalar result = 0;
    int j = gq_precision[precision].index;
    for (int i = 0; i < gq_precision[precision].num_points; i++) {
        f1->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, &q);
        result += gq_triangle[j].W * integrate_point_to_face(f, &q, f2, precision);
        j++;
    }
    return (result * f1->area);
}

