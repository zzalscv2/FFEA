#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include <math.h>

#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SecondOrderFunctions.h"
#include "SimulationParams.h"

class Face {
public:
    Face();

    ~Face();

    /* Pointers to the surface nodes that define this face */
    mesh_node *n[3];

    /* Pointer to the element this is a face of */
    tetra_element_linear *e;

    /* Van der Waals interaction type */
    int vdw_interaction_type;

    /* Initial, equilibrium area of this face */
    scalar area_0;

    /* Last calculated area (by most recent call to calc_area_normal_centroid()) */
    scalar area;

    /* Last calculated normal (by most recent call to calc_area_normal_centroid()) */
    vector3 normal;

    /* Last calculated centroid (by most recent call to calc_area_normal_centroid()) */
    vector3 centroid;

    /* Stores the current force applied to this face */
    vector3 force[3];


    /* Stores the natural (shape function) coords of the centroid of this face in the parent element */
    SecondOrderFunctions::stu centroid_stu;

    bool vdw_xz_interaction_flag;
    bool *vdw_bb_interaction_flag;

    /* vdw measurement info */
    int num_blobs;
    vector3 *vdw_bb_force;
    scalar *vdw_bb_energy;
    vector3 *vdw_xz_force;
    scalar vdw_xz_energy;

    Blob *daddy_blob;

    void init(tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, SecondOrderFunctions::stu centroid_stu, Blob *daddy_blob, SimulationParams *params);

    void init(mesh_node *n0, mesh_node *n1, mesh_node *n2, Blob *daddy_blob, SimulationParams *params);

    void set_vdw_interaction_type(int vdw_interaction_type);

    void calc_area_normal_centroid();

    /* Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3 */
    void barycentric_calc_point(scalar b1, scalar b2, scalar b3, vector3 *p);

    /* Returns the average electrostatic potential of this face */
    scalar average_phi();

    scalar get_normal_flux();

    void add_force_to_node(int i, vector3 *f);

    void add_force_to_node_atomic(int i, vector3 *f);

    void add_bb_vdw_force_to_record(vector3 *f, int other_blob_index);

    void add_bb_vdw_energy_to_record(scalar energy, int other_blob_index);

    void add_xz_vdw_force_to_record(vector3 *f);

    void add_xz_vdw_energy_to_record(scalar energy);

    void zero_force();

    void zero_vdw_bb_measurement_data();

    void zero_vdw_xz_measurement_data();

    void set_vdw_xz_interaction_flag(bool state);

    void set_vdw_bb_interaction_flag(bool state, int other_blob_index);

    bool is_vdw_active();

};

#endif
