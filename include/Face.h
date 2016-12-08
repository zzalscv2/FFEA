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

#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include <math.h>

#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SecondOrderFunctions.h"
#include "SimulationParams.h"
#include "CheckTetrahedraOverlap.h"
#include "VolumeIntersection.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"

class Face {
public:
    Face();

    ~Face();

    /** The first 3 pointers are to those surface nodes that define this face,
      * while the last one belongs to the opposite node of the linear tetrahedron **/
    mesh_node *n[4];

    /** Index of this face */
    int index;

    /** Pointer to the element this is a face of */
    tetra_element_linear *e;

    /** Van der Waals interaction type **/
    int vdw_interaction_type;

    /** Initial, equilibrium area of this face **/
    scalar area_0;

    /** Last calculated area (by most recent call to calc_area_normal_centroid()) **/
    scalar area;

    /** Last calculated normal (by most recent call to calc_area_normal_centroid()) **/
    vector3 normal;

    /** Last calculated centroid (by most recent call to calc_area_normal_centroid()) **/
    vector3 centroid;

    /** Stores the current force applied to this face.
      * The first 3 "vectors" are to those surface nodes that define this face,
      * while the last one belongs to the opposite node of the linear tetrahedron **/
    vector3 force[4];


    /** Stores the natural (shape function) coords of the centroid of this face in the parent element **/
    SecondOrderFunctions::stu centroid_stu;

    bool vdw_xz_interaction_flag;
    bool *vdw_bb_interaction_flag;
    bool kinetically_active;

    //@{
    /** vdw measurement info **/
    int num_blobs;
    vector3 *vdw_bb_force;
    scalar *vdw_bb_energy;
    vector3 *vdw_xz_force;
    scalar vdw_xz_energy;
    //@}

    /** Check whether the tetrahedron formed by this face an the opposite
      *   linear node does intersect with the corresponding tetrahedron in f2
      * Returns true if there is intersection.
      * Uses the "Fast Tetrahedron-Tetrahedron Overlap Algorithm".
      * It calls some private functions.
      **/
    bool checkTetraIntersection(Face *f2);


    /** Check whether the tetrahedron formed by this face an the opposite
      *   linear node does intersect with the corresponding tetrahedron in f2
      * Returns true if there is intersection.
      * Uses the "Fast Tetrahedron-Tetrahedron Overlap Algorithm".
      * It calls some private functions.
      * Altered to act periodically around box boundaries.
      **/
    bool checkTetraIntersection(Face *f2, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Get the volume that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      **/
    scalar getTetraIntersectionVolume(Face *f2);

    /** Check whether the tetrahedron formed by this face an the opposite
      *   linear node does intersect with the corresponding tetrahedron in f2.
      * If so, return the overlapping volume, otherwise return 0.
      * Uses the "Fast Tetrahedron-Tetrahedron Overlap Algorithm" for checking if 
      *  interaction occurs. 
      **/
    scalar checkTetraIntersectionAndGetVolume(Face *f2);


        /** Get the volume that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      **/
    scalar getTetraIntersectionVolume(Face *f2, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Get the volume and area that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      **/
    void getTetraIntersectionVolumeAndArea(Face *f2, geoscalar &vol, geoscalar &area);

        /** Get the volume and area that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      * Altered to act periodically around box boundaries.
      **/
    void getTetraIntersectionVolumeAndArea(Face *f2, geoscalar &vol, geoscalar &area, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);


    /** Get the volume that enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      *   In addition, return the gradient of this volume,
      *     with respect to the unit vector r.
      * It calls twice volumeIntersection, at volumeIntersection.h
      **/
    void getTetraIntersectionVolumeAndGradient(Face *f2, grr3 &r, geoscalar &vol, geoscalar &dVdr);

     /** Get the volume that enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      *   In addition, return the gradient of this volume,
      *     with respect to the unit vector r.
      * It calls volumeIntersection, at volumeIntersection.h
    * Altered to act periodically around box boundaries.
      **/
    void getTetraIntersectionVolumeAndGradient(Face *f2, grr3 &r, geoscalar &vol, geoscalar &dVdr, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

     /** Get the volume that enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      *   In addition, return the gradient of this volume,
      *     with respect to the unit vector r,
      *     and the two action points where the force is applied.
      * It calls twice volumeIntersection and lineFaceIntersectionPoint.
      **/
    bool getTetraIntersectionVolumeGradientAndShapeFunctions(Face *f2, grr3 (&r), geoscalar &vol, geoscalar &dVdr, grr4 (&phi1), grr4 (&phi2));

    Blob *daddy_blob;

    void init(int index, tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *oposite, SecondOrderFunctions::stu centroid_stu, Blob *daddy_blob, SimulationParams *params);

    void init(int index, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, Blob *daddy_blob, SimulationParams *params);

    void set_vdw_interaction_type(int vdw_interaction_type);

    void build_opposite_node();

    void calc_area_normal_centroid();
    void set_kinetic_state(bool state);

    vector3 * get_centroid();
    void print_centroid();
    void print_nodes();
    scalar get_area();

    /** Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3. Altered to act periodically around box boundaries. **/
    void barycentric_calc_point(scalar b1, scalar b2, scalar b3, vector3 *p);
    void barycentric_calc_point_f2(scalar b1, scalar b2, scalar b3, vector3 *p,scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Returns the average electrostatic potential of this face **/
    scalar average_phi();

    scalar get_normal_flux();

    void add_force_to_node(int i, vector3 *f);
    void add_force_to_node(int i, arr3 (&f));

    void add_force_to_node_atomic(int i, vector3 *f);

    void add_bb_vdw_force_to_record(vector3 *f, int other_blob_index);
    void add_bb_vdw_force_to_record(arr3 &f, int other_blob_index);

    void add_bb_vdw_energy_to_record(scalar energy, int other_blob_index);

    void add_xz_vdw_force_to_record(vector3 *f);

    void add_xz_vdw_energy_to_record(scalar energy);

    void zero_force();

    void zero_vdw_bb_measurement_data();

    void zero_vdw_xz_measurement_data();

    void set_vdw_xz_interaction_flag(bool state);

    void set_vdw_bb_interaction_flag(bool state, int other_blob_index);

    void vec3Vec3SubsToArr3Mod(Face *f2, arr3 (&w), scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    bool is_vdw_active();
    bool is_kinetic_active();

private:
    int stuff;

};

#endif
