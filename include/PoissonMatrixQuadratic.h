#ifndef POISSONMATRIXQUADRATIC_H_INCLUDED
#define POISSONMATRIXQUADRATIC_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_10X10 55

#define GRAD_PSI_1 0
#define GRAD_PSI_2 1
#define GRAD_PSI_3 2
#define GRAD_PSI_4 3
#define GRAD_PSI_5 4
#define GRAD_PSI_6 5
#define GRAD_PSI_7 6
#define GRAD_PSI_8 7
#define GRAD_PSI_9 8
#define GRAD_PSI_10 9

#define NUM_TET_GAUSS_QUAD_POINTS 14

#include "SecondOrderFunctions.h"

class PoissonMatrixQuadratic {
public:
    PoissonMatrixQuadratic();

    scalar * get_K_alpha_mem_loc(int i, int j);

    void build(mesh_node *n[10], scalar epsilon);

    scalar get_K_alpha_value(int i, int j);

private:
    scalar K_alpha[NUM_ELEMENTS_LOWER_TRIANGULAR_10X10];

    struct tetrahedron_gauss_point {
        scalar W;
        scalar eta[4];
    };

    void add_grad_dot_products(vector3 grad_psi[10], scalar det_J, scalar weight);

    scalar grad_dot(vector3 *grad_psi_i, vector3 *grad_psi_j);

    void zero();

    void scale(scalar factor);
};


#endif
