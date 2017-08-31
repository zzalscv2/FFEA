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

#include "VdW_solver.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

// const scalar VdW_solver::phi_f[4] = { 0.25, 0.25, 0.25, 0.25};

const int VdW_solver::adjacent_cell_lookup_table[27][3] = {
    {-1, -1, -1},
    {-1, -1, 0},
    {-1, -1, +1},
    {-1, 0, -1},
    {-1, 0, 0},
    {-1, 0, +1},
    {-1, +1, -1},
    {-1, +1, 0},
    {-1, +1, +1},
    {0, -1, -1},
    {0, -1, 0},
    {0, -1, +1},
    {0, 0, -1},
    {0, 0, 0},
    {0, 0, +1},
    {0, +1, -1},
    {0, +1, 0},
    {0, +1, +1},
    {+1, -1, -1},
    {+1, -1, 0},
    {+1, -1, +1},
    {+1, 0, -1},
    {+1, 0, 0},
    {+1, 0, +1},
    {+1, +1, -1},
    {+1, +1, 0},
    {+1, +1, +1}
};

const VdW_solver::tri_gauss_point VdW_solver::gauss_points[] = {
    // Weight, eta1, eta2, eta3
    {   0.333333333333333,
        {0.666666666666667, 0.166666666666667, 0.166666666666667}
    },
    {   0.333333333333333,
        {0.166666666666667, 0.666666666666667, 0.166666666666667}
    },
    {   0.333333333333333,
        {0.166666666666667, 0.166666666666667, 0.666666666666667}
    }

    /*
                                            {0.109951743655322,     {0.816847572980459, 0.091576213509771, 0.091576213509771}},
                                            {0.109951743655322,     {0.091576213509771, 0.816847572980459, 0.091576213509771}},
                                            {0.109951743655322,     {0.091576213509771, 0.091576213509771, 0.816847572980459}},
                                            {0.223381589678011,     {0.108103018168070, 0.445948490915965, 0.445948490915965}},
                                            {0.223381589678011,     {0.445948490915965, 0.108103018168070, 0.445948490915965}},
                                            {0.223381589678011,     {0.445948490915965, 0.445948490915965, 0.108103018168070}}
     */

    /*
                                            {0.050844906370207,     {0.873821971016996, 0.063089014491502, 0.063089014491502}},
                                            {0.050844906370207,     {0.063089014491502, 0.873821971016996, 0.063089014491502}},
                                            {0.050844906370207,     {0.063089014491502, 0.063089014491502, 0.873821971016996}},
                                            {0.116786275726379,     {0.501426509658179, 0.249286745170910, 0.249286745170910}},
                                            {0.116786275726379,     {0.249286745170910, 0.501426509658179, 0.249286745170910}},
                                            {0.116786275726379,     {0.249286745170910, 0.249286745170910, 0.501426509658179}},
                                            {0.082851075618374,     {0.636502499121399, 0.310352451033785, 0.053145049844816}},
                                            {0.082851075618374,     {0.310352451033785, 0.053145049844816, 0.636502499121399}},
                                            {0.082851075618374,     {0.053145049844816, 0.636502499121399, 0.310352451033785}},
                                            {0.082851075618374,     {0.636502499121399, 0.053145049844816, 0.310352451033785}},
                                            {0.082851075618374,     {0.310352451033785, 0.636502499121399, 0.053145049844816}},
                                            {0.082851075618374,     {0.053145049844816, 0.310352451033785, 0.636502499121399}}
     */
};

VdW_solver::VdW_solver() {
    total_num_surface_faces = 0;
    surface_face_lookup = NULL;
    box_size.x = 0;
    box_size.y = 0;
    box_size.z = 0;
    num_blobs = 0;
    fieldenergy = NULL;
    vdw_type = VDW_TYPE_UNDEFINED;
}

VdW_solver::~VdW_solver() {
    total_num_surface_faces = 0;
    surface_face_lookup = NULL;
    box_size.x = 0;
    box_size.y = 0;
    box_size.z = 0;
    for(int i = 0; i < num_blobs; ++i) {
        delete[] fieldenergy[i];
    }
    delete[] fieldenergy;
    fieldenergy = NULL;
    num_blobs = 0;
    vdw_type = VDW_TYPE_UNDEFINED;
}

int VdW_solver::init(NearestNeighbourLinkedListCube *surface_face_lookup, vector3 *box_size, LJ_matrix *lj_matrix, scalar &vdw_steric_factor, int num_blobs, int inc_self_vdw, string vdw_type_string, scalar &vdw_steric_dr, int calc_kinetics, bool working_w_static_blobs) {
    this->surface_face_lookup = surface_face_lookup;
    this->box_size.x = box_size->x;
    this->box_size.y = box_size->y;
    this->box_size.z = box_size->z;

    this->lj_matrix = lj_matrix;

    this->inc_self_vdw = inc_self_vdw;
    this->steric_factor = vdw_steric_factor;
    this->steric_dr = vdw_steric_dr;
    if (vdw_type_string == "lennard-jones")
        vdw_type = VDW_TYPE_LJ;
    else if (vdw_type_string == "steric")
        vdw_type = VDW_TYPE_STERIC;
    else if (vdw_type_string == "ljsteric")
        vdw_type = VDW_TYPE_LJSTERIC;


    // And some measurement stuff it should know about
    this->num_blobs = num_blobs;
    this->calc_kinetics = calc_kinetics;
    this->working_w_static_blobs = working_w_static_blobs;
    fieldenergy = new scalar*[num_blobs];
    if (fieldenergy == NULL) FFEA_ERROR_MESSG("Failed to allocate fieldenergy in VdW_solver::init\n");
    for(int i = 0; i < num_blobs; ++i) {
        fieldenergy[i] = new scalar[num_blobs];
        if (fieldenergy[i] == NULL) FFEA_ERROR_MESSG("Failed to allocate fieldenergy[%d] in VdW_solver::init\n", i);
    }
    return FFEA_OK;
}

/**  Zero measurement stuff, AKA fieldenergy */
void VdW_solver::reset_fieldenergy() {
    for(int i = 0; i < num_blobs; ++i) {
        for(int j = 0; j < num_blobs; ++j) {
            fieldenergy[i][j] = 0.0;
        }
    }
}

/** Solve VdW */ 
int VdW_solver::solve(scalar *blob_corr) {

    LinkedListNode<Face> *l_i = NULL;
    LinkedListNode<Face> *l_j = NULL;
    Face *f_i, *f_j;
    int c;
    total_num_surface_faces = surface_face_lookup->get_pool_size();
    //total_num_surface_faces = surface_face_lookup->get_stack_size();

    reset_fieldenergy();
    int motion_state_i;

    /* For each face, calculate the interaction with all other relevant faces and add the contribution to the force on each node, storing the energy contribution to "blob-blob" (bb) interaction energy.*/
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(blob_corr) private(c, l_i, l_j, f_i, f_j, motion_state_i) schedule(dynamic, 1) // OMP-GHL
#endif
    for (int i = 0; i < total_num_surface_faces; i++) {

        // get the ith face
        l_i = surface_face_lookup->get_from_pool(i);
        if ((calc_kinetics == 1) && (!l_i->obj->is_kinetic_active()) ) {
            continue;
        }
        f_i = l_i->obj;
        if (working_w_static_blobs) motion_state_i = f_i->daddy_blob->get_motion_state();
        int l_index_i = l_i->index;

        // Calculate this face's interaction with all faces in its cell and the 26 adjacent cells (3^3 = 27 cells)
        // Remember to check that the face is not interacting with itself or connected faces
        for (c = 0; c < 27; c++) {
            l_j = surface_face_lookup->get_top_of_stack(
                      l_i->x + adjacent_cell_lookup_table[c][0],
                      l_i->y + adjacent_cell_lookup_table[c][1],
                      l_i->z + adjacent_cell_lookup_table[c][2]);
            while (l_j != NULL) {
                if (consider_interaction(f_i, l_index_i, motion_state_i, l_j, blob_corr)) {
                    do_interaction(f_i, l_j->obj, blob_corr);
                }
                l_j = l_j->next;
            }
        }
    }
   // exit(0);
    return FFEA_OK;
}


/* Allow protein VdW interactions along the top and bottom x-z planes */
int VdW_solver::solve_sticky_wall(scalar h) {
    int Nx = 0, Ny = 0, Nz = 0;
    surface_face_lookup->get_dim(&Nx, &Ny, &Nz);
    LinkedListNode<Face> *l_j = NULL;
    Face *f_j = NULL;
    for (int y = 0; y < Ny; y += Ny - 1) {
        for (int z = 0; z < Nz; z++) {
            for (int x = 0; x < Nx; x++) {
                l_j = surface_face_lookup->get_top_of_stack(x, y, z);
                while (l_j != NULL) {
                    f_j = l_j->obj;
                    f_j->set_vdw_xz_interaction_flag(true);
                    do_sticky_xz_interaction(f_j, (y == 0), h * Ny);
                    l_j = l_j->next;
                }
            }
        }
    }
    return FFEA_OK;
}

void VdW_solver::do_lj_interaction(Face *f1, Face *f2, scalar *blob_corr) {

    int f1_daddy_blob_index = f1->daddy_blob->blob_index;
    int f2_daddy_blob_index = f2->daddy_blob->blob_index;

    //printf("Centroid 1: %f %f %f\n", f1->centroid.x * (mesoDimensions::length / 1e-9), f1->centroid.y * (mesoDimensions::length / 1e-9), f1->centroid.z * (mesoDimensions::length / 1e-9));
    //printf("Centroid 2: %f %f %f\n", f2->centroid.x * (mesoDimensions::length / 1e-9), f2->centroid.y * (mesoDimensions::length / 1e-9), f2->centroid.z * (mesoDimensions::length / 1e-9));
    //printf("Separation: %f %f %f\n", (f1->centroid.x - f2->centroid.x) * (mesoDimensions::length / 1e-9), (f1->centroid.y - f2->centroid.y) * (mesoDimensions::length / 1e-9), (f1->centroid.z - f2->centroid.z) * (mesoDimensions::length / 1e-9));

    // Get the interaction LJ parameters for these two face types
    scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
    lj_matrix->get_LJ_params(f1->vdw_interaction_type, f2->vdw_interaction_type, &vdw_eps, &vdw_r_eq);
    //printf("VdW Params: %f nm %e \n", vdw_r_eq * (mesoDimensions::length / 1e-9), vdw_eps * (mesoDimensions::Energy / mesoDimensions::area) * (1.0/ mesoDimensions::area));
    vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_quad_points];
    vector3 force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

    // Convert all area coordinate gauss points to cartesian
    if (blob_corr == NULL) {
        for (int i = 0; i < num_tri_gauss_quad_points; i++) {
            f1->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
            f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
        }
    } else {
        for (int i = 0; i < num_tri_gauss_quad_points; i++) {
            f1->barycentric_calc_point_f2(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i],blob_corr, f1_daddy_blob_index, f2_daddy_blob_index);
            f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
        }
    } 

    // Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
    // Also calculate energy whilst looping through face points
    scalar energy = 0.0;

    if (vdw_type == VDW_TYPE_LJSTERIC) calc_ljinterpolated_force_pair_matrix(force_pair_matrix, 
             p, q, vdw_r_eq, vdw_eps, energy);
    else if (vdw_type == VDW_TYPE_LJ) calc_lj_force_pair_matrix(force_pair_matrix, 
             p, q, vdw_r_eq, vdw_eps, energy);

    scalar ApAq = f1->area * f2->area;
    energy *= ApAq;

    // Store the measurement
	/*for (int k = 0; k < num_tri_gauss_quad_points; k++) {
		for (int l = 0; l < num_tri_gauss_quad_points; l++) {
			printf("%e  ", force_pair_matrix[k][l].z);
		}
		printf("\n");
	}
	exit(0);*/
    #pragma omp critical
    {
        fieldenergy[f1_daddy_blob_index][f2_daddy_blob_index] += energy;
        for (int j = 0; j < 3; j++) {
            vector3 force1, force2;
            force1.assign( 0, 0, 0 );
            force2.assign( 0, 0, 0 );
            for (int k = 0; k < num_tri_gauss_quad_points; k++) {
                for (int l = 0; l < num_tri_gauss_quad_points; l++) {
                    scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
                    scalar d = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];
                    //printf("c = %e, %e, %e, %e\n", c, gauss_points[k].W, gauss_points[l].W, gauss_points[l].eta[j]);
                    force1.x += c * force_pair_matrix[k][l].x;
                    force1.y += c * force_pair_matrix[k][l].y;
                    force1.z += c * force_pair_matrix[k][l].z;

                    force2.x -= d * force_pair_matrix[l][k].x;
                    force2.y -= d * force_pair_matrix[l][k].y;
                    force2.z -= d * force_pair_matrix[l][k].z;
                }
            }
            force1.x *= ApAq;
            force1.y *= ApAq;
            force1.z *= ApAq;
            f1->add_force_to_node(j, &force1);
            // f1->add_bb_vdw_force_to_record(&force1, f2->daddy_blob->blob_index); // DEPRECATED
            			//	printf("1:: %d %e %e %e\n", f1->index, force1.x, force1.y, force1.z);
					//printf("2:: %d %e %e %e\n", f2->index, force2.x, force2.y, force2.z);

            force2.x *= ApAq;
            force2.y *= ApAq;
            force2.z *= ApAq;
            f2->add_force_to_node(j, &force2);
            // f2->add_bb_vdw_force_to_record(&force2, f1->daddy_blob->blob_index); // DEPRECATED
        } // end updating face nodes.
    } // end of critical


}

/**Alters interaction calculations to apply periodic boundary conditions*/
void VdW_solver::do_interaction(Face *f1, Face *f2, scalar *blob_corr) {

    do_lj_interaction(f1, f2, blob_corr);

}

bool VdW_solver::consider_interaction(Face *f_i, int l_index_i, int motion_state_i, LinkedListNode<Face> *l_j, scalar *blob_corr) {

    bool interaction_needed = false;
    if (l_index_i < l_j->index) {
        if ((inc_self_vdw == 1) or ( (inc_self_vdw == 0 ) and (f_i->daddy_blob != l_j->obj->daddy_blob))) {

            if((working_w_static_blobs == false) || (motion_state_i == FFEA_BLOB_IS_DYNAMIC or l_j->obj->daddy_blob->get_motion_state() == FFEA_BLOB_IS_DYNAMIC)) {
                interaction_needed = true;
            }
        }
    }

    if (interaction_needed) {

        Face *f_j = l_j->obj;
        // 1 - Check that faces are facing each other, if not then they are not interacting
	//cout << f_i->normal[0]*f_j->normal[0] + f_i->normal[1]*f_j->normal[1] + f_i->normal[2]*f_j->normal[2] << endl;
        if ( (f_i->normal[0]*f_j->normal[0] +
                f_i->normal[1]*f_j->normal[1] +
                f_i->normal[2]*f_j->normal[2]) > ffea_const::zero ) return false;


        if ( vdw_type != VDW_TYPE_LJ ) {
            // 2 - Two more checks:
            // 2.1 - Check that faces are in front of each other
            //     - Robin suspected this was leading to unstabilities for the LJ case.
            vector3 sep;
            if (blob_corr == NULL) {
                sep.assign ( f_j->centroid.x - f_i->centroid.x, f_j->centroid.y - f_i->centroid.y, f_j->centroid.z - f_i->centroid.z );
            } else {
                sep.assign ( f_j->centroid.x - f_i->centroid.x-blob_corr[f_i->daddy_blob->blob_index*(num_blobs)*3 + f_j->daddy_blob->blob_index*3],f_j->centroid.y - f_i->centroid.y-blob_corr[f_i->daddy_blob->blob_index*(num_blobs)*3 + f_j->daddy_blob->blob_index*3+1],f_j->centroid.z - f_i->centroid.z-blob_corr[f_i->daddy_blob->blob_index*(num_blobs)*3 + f_j->daddy_blob->blob_index*3+2] );
            }
            if ((arr3arr3DotProduct<scalar,arr3>(sep.data, f_i->normal.data) < 0) &&
                    arr3arr3DotProduct<scalar,arr3>(sep.data, f_j->normal.data) > 0) return false;


            // 2.2 - Check that no nodes are shared,
            //     only in the case that faces belong to the same blob:
            if ((inc_self_vdw == 1) && (f_i->daddy_blob == f_j->daddy_blob)) {
                if (f_i->n[3] == f_j->n[3]) {
                    return false;
                }
                for (int i=0; i<4; i++) {
                    int in_i = f_i->n[i]->index;
                    for (int j=0; j<4; j++) {
                        if (f_j->n[j]->index == in_i) {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return interaction_needed;

}

void VdW_solver::do_sticky_xz_interaction(Face *f, bool bottom_wall, scalar dim_y) {

    scalar y_wall = 0; //-vdw_r_eq;
    if (bottom_wall == false) {
        y_wall = dim_y; // + vdw_r_eq;
    }

    // Check that face is facing wall, if not then it should not be interacting
    if ((bottom_wall == true && f->normal.y > 0) || (bottom_wall == false && f->normal.y < 0)) {
        return;
    }

    // Get the interaction LJ parameters for these two face types
    scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
    lj_matrix->get_LJ_params(f->vdw_interaction_type, f->vdw_interaction_type, &vdw_eps, &vdw_r_eq);


    vector3 p[num_tri_gauss_quad_points];
    scalar force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

    // Convert all area coordinate gauss points to cartesian
    for (int i = 0; i < num_tri_gauss_quad_points; i++) {
        f->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
    }

    // Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
    // Also calculate the energy whilst looping through face points
    scalar energy = 0.0;
    for (int k = 0; k < num_tri_gauss_quad_points; k++) {
        for (int l = k; l < num_tri_gauss_quad_points; l++) {
            scalar mag_r = p[k].y - y_wall;

            scalar force_mag = 12 * pow(vdw_r_eq, 6) * vdw_eps * (pow(mag_r, -7) - pow(vdw_r_eq, 6) * pow(mag_r, -13));
            energy += pow(vdw_r_eq, 6) * vdw_eps * (pow(vdw_r_eq, 6) * pow(mag_r, -12) - 2 * pow(mag_r, -6));
            force_mag *= -1;

            force_pair_matrix[k][l] = force_mag;
            force_pair_matrix[l][k] = force_pair_matrix[k][l];
        }
    }

    // Record energy with xz plane
    scalar Asq = f->area * f->area;
    energy *= Asq;
    // f->add_xz_vdw_energy_to_record(energy); // DEPRECATED

    for (int j = 0; j < 3; j++) {
        vector3 force;
        force.assign( 0, 0, 0 );
        for (int k = 0; k < num_tri_gauss_quad_points; k++) {
            for (int l = 0; l < num_tri_gauss_quad_points; l++) {
                scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
                force.y += c * force_pair_matrix[k][l];
            }
        }
        force.y *= Asq;
        f->add_force_to_node(j, &force);
        // f->add_xz_vdw_force_to_record(&force); // DEPRECATED
    }

}

bool VdW_solver::do_steric_interaction(Face *f1, Face *f2, scalar *blob_corr) {

    int f1_daddy_blob_index = f1->daddy_blob->blob_index;
    int f2_daddy_blob_index = f2->daddy_blob->blob_index;

    // //  Working version for F = k*dV/dr // //
    geoscalar vol;
    grr3 dVdr;
    grr4 phi1, phi2;

    if (blob_corr == NULL) {
      if (!f1->getTetraIntersectionVolumeTotalGradientAndShapeFunctions(f2, steric_dr, dVdr, vol, phi1, phi2)) return false;
    } else {
      if (!f1->getTetraIntersectionVolumeTotalGradientAndShapeFunctions(f2, steric_dr, dVdr, vol, phi1, phi2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index)) return false;
    }

    vol *= steric_factor;

    // Force is proportional to the gradient, i. e.:
    arr3Resize<geoscalar,grr3>(steric_factor, dVdr);

    grr3 ftmp1, ftmp2;
    #pragma omp critical
    {
        // Store the measurement
        fieldenergy[f1_daddy_blob_index][f2_daddy_blob_index] += vol;
        // Finally, apply the force onto the nodes:
        for (int j = 0; j < 4; j++) {
            arr3Resize2<geoscalar,grr3>(phi1[j], dVdr, ftmp1);
            f1->add_force_to_node(j, ftmp1);
            // f1->add_bb_vdw_force_to_record(ftmp1, f2->daddy_blob->blob_index); // DEPRECATED

            arr3Resize2<geoscalar,grr3>(ffea_const::mOne*phi2[j], dVdr, ftmp2);
            f2->add_force_to_node(j, ftmp2);
            // f2->add_bb_vdw_force_to_record(ftmp2, f1->daddy_blob->blob_index); // DEPRECATED
        }
    }

    /* // //  Working version for F = k // //
    geoscalar vol, dVdr;
    grr3 force1, force2; //, n1_b;
    //  Then, check whether the tetrahedra intersect,
    //    and if so, get the volume:
    scalar vol = f1->checkTetraIntersectionAndGetVolume(f2);
    if ( vol < ffea_const::threeErr ) return;

    // Choose the force line
    // as the line passing through the elements CMs.
    arr3 force1, force2, cm1, cm2; //, n1_b;
    arr3Initialise<arr3>(cm1);
    arr3Initialise<arr3>(cm2);
    for (int i=0; i<4; i++) {
      cm1[0] += f1->n[i]->pos.x;
      cm1[1] += f1->n[i]->pos.y;
      cm1[2] += f1->n[i]->pos.z;
      cm2[0] += f2->n[i]->pos.x;
      cm2[1] += f2->n[i]->pos.y;
      cm2[2] += f2->n[i]->pos.z;
    }
    arr3Resize<scalar,arr3>(0.25,cm1);
    arr3Resize<scalar,arr3>(0.25,cm2);
    arr3arr3Substract<scalar,arr3>(cm2, cm1, force2);
    //printf("**********\n Blob %d to Blob %d\n face %d to face %d\ndist in x is %f\ndist in y is %f\ndist in z is %f\n",f1->daddy_blob->blob_index,f2->daddy_blob->blob_index,f1->index, f2->index,force2[0],force2[1],force2[2]);
    arr3Normalise<scalar,arr3>(force2); // that is the direction of the force for f2 (backwards).

    // Store the measurement
    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;

    // Force is proportional to the gradient, i. e.:
    arr3Resize<scalar,arr3>(steric_factor, force2);
    arr3Resize2<scalar,arr3>(ffea_const::mOne, force2, force1);

    // Finally, apply the force onto the nodes:
    for (int j = 0; j < 4; j++) {
      f1->add_force_to_node(j, force1);
      f1->add_bb_vdw_force_to_record(force1, f2->daddy_blob->blob_index);
      f2->add_force_to_node(j, force2);
      f2->add_bb_vdw_force_to_record(force2, f1->daddy_blob->blob_index);
    } */

    return true;

}


scalar VdW_solver::distance2(vector3 &p, vector3 &q) {
    scalar dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
    return dx * dx + dy * dy + dz * dz;
}

scalar VdW_solver::dot(vector3 *p, vector3 *q) {
    return p->x * q->x + p->y * q->y + p->z * q->z;
}

scalar VdW_solver::dot_with_normal(vector3 *p, vector3 *q, vector3 *n) {
    return (q->x - p->x) * n->x + (q->y - p->y) * n->y + (q->z - p->z) * n->z;
}

scalar VdW_solver::minimum_image(scalar delta, scalar size) {
    if (fabs(delta) > size * .5) {
        return size - delta;
    }

    return delta;
}

scalar VdW_solver::get_field_energy(int index0, int index1) {

    // Sum over all field
    if(index0 == -1 || index1 == -1) {
        scalar energy = 0.0;
        for(int i = 0; i < num_blobs; ++i) {
            for(int j = 0; j < num_blobs; ++j) {
                energy += fieldenergy[i][j];
            }
        }

        return energy;

    } else if (index0 == index1) {
        return fieldenergy[index0][index1];
    } else {

        // Order of blob indices is unknown in the calculations, so must add
        return fieldenergy[index0][index1] + fieldenergy[index1][index0];
    }
}

void VdW_solver::calc_lj_force_pair_matrix(vector3 (&force_pair_matrix)[num_tri_gauss_quad_points][num_tri_gauss_quad_points],
        vector3 (&p)[num_tri_gauss_quad_points], vector3 (&q)[num_tri_gauss_quad_points],
        scalar &vdw_r_eq, scalar &vdw_eps, scalar &energy) {

    scalar mag_r, force_mag, e;

    scalar vdw_r_eq_2 = vdw_r_eq * vdw_r_eq;
    scalar vdw_r_eq_4 = vdw_r_eq_2 * vdw_r_eq_2;
    scalar vdw_r_eq_6 = vdw_r_eq_4 * vdw_r_eq_2;
   
    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
        mag_r = sqrt(distance2(p[k], q[k])); 
        calc_lj_factors(mag_r, k, k, vdw_eps, vdw_r_eq_6, force_mag, e); 
        energy += e; 
        force_pair_matrix[k][k].x = force_mag * ((p[k].x - q[k].x) / mag_r);
        force_pair_matrix[k][k].y = force_mag * ((p[k].y - q[k].y) / mag_r);
        force_pair_matrix[k][k].z = force_mag * ((p[k].z - q[k].z) / mag_r);

        for(int l = k+1; l < num_tri_gauss_quad_points; l++) {
            mag_r = sqrt(distance2(p[k], q[l])); 
            calc_lj_factors(mag_r, k, l, vdw_eps, vdw_r_eq_6, force_mag, e); 

            energy += 2*e;

            force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
            force_pair_matrix[k][l].y = force_mag * ((p[k].y - q[l].y) / mag_r);
            force_pair_matrix[k][l].z = force_mag * ((p[k].z - q[l].z) / mag_r);

            force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
            force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
            force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
        }
    }
}

/** Given (mag_r), get LJ force magnitude (force_mag) and energy (e) */
void VdW_solver::calc_lj_factors(scalar &mag_r, int index_k, int index_l, scalar &vdw_eps, scalar &vdw_r_eq_6,
                                 scalar &force_mag, scalar &e) {

    scalar mag_ri,  mag_ri_2, mag_ri_4, mag_ri_6, mag_ri_7;
    scalar vdw_fac_6;

    mag_ri = 1./mag_r;
    mag_ri_2 = mag_ri * mag_ri;
    mag_ri_4 = mag_ri_2 * mag_ri_2;
    mag_ri_6 = mag_ri_4 * mag_ri_2;
    mag_ri_7 = mag_ri_6 * mag_ri;
    vdw_fac_6 = vdw_r_eq_6 * mag_ri_6;
    force_mag = 12 * mag_ri_7 * vdw_r_eq_6 * vdw_eps * (vdw_fac_6 - 1);
    e = gauss_points[index_k].W * gauss_points[index_l].W *
                   vdw_eps * vdw_fac_6 * (vdw_fac_6 - 2 );

}

/** Given (mag_r), get LJ_interpolated force magnitude (force_mag) and energy (e) */
void VdW_solver::calc_ljinterpolated_factors(scalar &mag_r, int index_k, int index_l, scalar &vdw_eps, scalar &vdw_r_eqi,
                                 scalar &force_mag, scalar &e) {

    scalar vdw_fac = mag_r * vdw_r_eqi;
    scalar vdw_fac_2 = vdw_fac * vdw_fac;

    e = gauss_points[index_k].W * gauss_points[index_l].W *
                    vdw_eps * vdw_fac_2 * (2 * vdw_fac - 3);
    force_mag = 6 * vdw_eps * vdw_r_eqi * vdw_fac * (1 - vdw_fac);

}

void VdW_solver::calc_ljinterpolated_force_pair_matrix(vector3 (&force_pair_matrix)[num_tri_gauss_quad_points][num_tri_gauss_quad_points],
        vector3 (&p)[num_tri_gauss_quad_points], vector3 (&q)[num_tri_gauss_quad_points],
        scalar &vdw_r_eq, scalar &vdw_eps, scalar &energy) {

    scalar mag_r, e, force_mag;

    scalar vdw_r_eq_2 = vdw_r_eq * vdw_r_eq;
    scalar vdw_r_eq_4 = vdw_r_eq_2 * vdw_r_eq_2;
    scalar vdw_r_eq_6 = vdw_r_eq_4 * vdw_r_eq_2;

    scalar vdw_r_eqi = 1.0 / vdw_r_eq;

    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
        mag_r = sqrt(distance2(p[k], q[k]));
        if(mag_r < vdw_r_eq) 
           calc_ljinterpolated_factors(mag_r, k, k, vdw_eps, vdw_r_eqi, force_mag, e);
        else 
           calc_lj_factors(mag_r, k, k, vdw_eps, vdw_r_eq_6, force_mag, e);
        energy += e;
        force_pair_matrix[k][k].x = force_mag * ((p[k].x - q[k].x) / mag_r);
        force_pair_matrix[k][k].y = force_mag * ((p[k].y - q[k].y) / mag_r);
        force_pair_matrix[k][k].z = force_mag * ((p[k].z - q[k].z) / mag_r);

        for(int l = k+1; l < num_tri_gauss_quad_points; l++) {

            mag_r = sqrt(distance2(p[k], q[l]));

            if(mag_r < vdw_r_eq) 
                calc_ljinterpolated_factors(mag_r, k, l, vdw_eps, vdw_r_eqi, force_mag, e);
            else
                calc_lj_factors(mag_r, k, l, vdw_eps, vdw_r_eq_6, force_mag, e);

            energy += 2*e;

            force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
            force_pair_matrix[k][l].y = force_mag * ((p[k].y - q[l].y) / mag_r);
            force_pair_matrix[k][l].z = force_mag * ((p[k].z - q[l].z) / mag_r);

            force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
            force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
            force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
        }
    }

}
