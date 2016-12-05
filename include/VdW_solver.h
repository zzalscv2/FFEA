#ifndef VDW_SOLVER_H_INCLUDED
#define VDW_SOLVER_H_INCLUDED

#include <math.h>

#include "FFEA_return_codes.h"
#include "NearestNeighbourLinkedListCube.h"
#include "LJ_matrix.h"
#include "Blob.h"

class VdW_solver {
public:
    VdW_solver();

    ~VdW_solver();

    int init(NearestNeighbourLinkedListCube *surface_face_lookup, vector3 *box_size, LJ_matrix *lj_matrix, scalar &vdw_steric_factor, int num_blobs, int inc_self_vdw);

    int solve();

    int solve(scalar * blob_corr);

    /** Allow protein VdW interactions along the top and bottom x-z planes */
    int solve_sticky_wall(scalar h);

    scalar get_field_energy(int i, int j);

protected:
    int total_num_surface_faces;
    NearestNeighbourLinkedListCube *surface_face_lookup;

    vector3 box_size;
    LJ_matrix *lj_matrix;

    scalar **fieldenergy;
    int num_blobs;
    int inc_self_vdw;  ///< whether to include interactions between faces within the same blob, or not.
    struct adjacent_cell_lookup_table_entry {
        int ix, iy, iz;
    };

    struct tri_gauss_point {
        scalar W;
        scalar eta[3];
    };

    virtual void do_interaction(Face *f1, Face *f2);

    virtual void do_interaction(Face *f1, Face *f2, scalar * blob_corr);

    scalar steric_factor; ///< Proportionality factor to the Steric repulsion.
    // static const scalar phi_f[4]; ///< shape function for the center of the "element"

    void do_sticky_xz_interaction(Face *f, bool bottom_wall, scalar dim_y);

    scalar distance2(vector3 *p, vector3 *q);

    scalar dot(vector3 *p, vector3 *q);

    scalar dot_with_normal(vector3 *p, vector3 *q, vector3 *n);

    scalar minimum_image(scalar delta, scalar size);
};

#endif
