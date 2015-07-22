#ifndef MASSMATRIXQUADRATIC_H_INCLUDED
#define MASSMATRIXQUADRATIC_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_10X10 55

#define NUM_SHAPE_FUNCTIONS 10

#define NUM_TET_GAUSS_QUAD_POINTS 14

#include "SecondOrderFunctions.h"

class MassMatrixQuadratic {
public:
    MassMatrixQuadratic();

    scalar * get_M_alpha_mem_loc(int i, int j);

    void build(mesh_node *n[10]);

    scalar get_M_alpha_value(int i, int j);

private:
    scalar M_alpha[NUM_ELEMENTS_LOWER_TRIANGULAR_10X10];

    struct tetrahedron_gauss_point {
        scalar W;
        scalar eta[4];
    };

    void add_psi_dot_products(scalar psi[10], scalar det_J, scalar weight);

    void zero();
};


#endif
