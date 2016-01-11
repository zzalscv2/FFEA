#include "VdW_solver.h"

VdW_solver::VdW_solver() {
    total_num_surface_faces = 0;
    surface_face_lookup = NULL;
    box_size.x = 0;
    box_size.y = 0;
    box_size.z = 0;
}

VdW_solver::~VdW_solver() {
    total_num_surface_faces = 0;
    surface_face_lookup = NULL;
    box_size.x = 0;
    box_size.y = 0;
    box_size.z = 0;
}

int VdW_solver::init(NearestNeighbourLinkedListCube *surface_face_lookup, vector3 *box_size, LJ_matrix *lj_matrix) {
    this->surface_face_lookup = surface_face_lookup;
    this->box_size.x = box_size->x;
    this->box_size.y = box_size->y;
    this->box_size.z = box_size->z;

    this->lj_matrix = lj_matrix;
    return FFEA_OK;
}

int VdW_solver::solve() {
    const struct adjacent_cell_lookup_table_entry adjacent_cell_lookup_table[27] ={
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

    LinkedListNode<Face> *l_i = NULL;
    LinkedListNode<Face> *l_j = NULL;
    Face *f_i, *f_j;

    total_num_surface_faces = surface_face_lookup->get_pool_size();

    /* For each face, calculate the interaction with all other relevant faces and add the contribution to the force on each node, storing the energy contribution to "blob-blob" (bb) interaction energy.*/ 
#ifdef USE_OPENMP
#pragma omp parallel for private(l_i, l_j, f_i, f_j) 
#endif
    for (int i = 0; i < total_num_surface_faces; i++) {

        // get the ith face
        l_i = surface_face_lookup->get_from_pool(i);
        f_i = l_i->obj;

        // Calculate this face's interaction with all faces in its cell and the 26 adjacent cells (3^3 = 27 cells)
        // Remember to check that the face is not interacting with itself or connected faces
        for (int c = 0; c < 27; c++) {
            l_j = surface_face_lookup->get_top_of_stack( 
                    l_i->x + adjacent_cell_lookup_table[c].ix,
                    l_i->y + adjacent_cell_lookup_table[c].iy,
                    l_i->z + adjacent_cell_lookup_table[c].iz);
            while (l_j != NULL) {
                if (l_i->index != l_j->index) {
                    f_j = l_j->obj;
                    if (f_i->daddy_blob != f_j->daddy_blob) {
                        //printf("(%d %d)\n", l_i->index, l_j->index);
                        //if((l_i->index == 1 && l_j->index == 5) || (l_i->index == 5 && l_j->index == 1))
                        f_i->set_vdw_bb_interaction_flag(true, f_j->daddy_blob->blob_index);
                        f_j->set_vdw_bb_interaction_flag(true, f_i->daddy_blob->blob_index);
                        do_interaction(f_i, f_j);
                        // do_volumeExclusion(f_i, f_j);
                    }
                }
                l_j = l_j->next;
            }
        }

    }

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


void VdW_solver::do_interaction(Face *f1, Face *f2) {
    // First check two things (either of which results in not having to calculate anything):
    // Check that faces are facing each other, if not then they are not interacting
    if (dot(&f1->normal, &f2->normal) > 0) {
        //				printf("DENIED NORMAL\n");
        return;
    }

    /* Robin suggested that this could lead to unstabilities for bad meshes... */ 
    // Check that faces are in front of each other
    //			vector3 sep = {f2->centroid.x - f1->centroid.x, f2->centroid.y - f1->centroid.y, f2->centroid.z - f1->centroid.z};
    //			if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
    //				printf("DENIED BEHIND\n");
    //				return;
    //			}

    // If faces are more than a few nanometres apart, don't bother calculating the force (it will be tiny)
    //			if(distance2(&f1->centroid, &f2->centroid) > (vdw_r_eq * 6 * vdw_r_eq * 6)) {
    //				return;
    //			}

    // Get the interaction LJ parameters for these two face types
    scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
    lj_matrix->get_LJ_params(f1->vdw_interaction_type, f2->vdw_interaction_type, &vdw_eps, &vdw_r_eq);

    const int num_tri_gauss_quad_points = 3;

    const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] ={
        // Weight, eta1, eta2, eta3
        {0.333333333333333,
            {0.666666666666667, 0.166666666666667, 0.166666666666667}},
        {0.333333333333333,
            {0.166666666666667, 0.666666666666667, 0.166666666666667}},
        {0.333333333333333,
            {0.166666666666667, 0.166666666666667, 0.666666666666667}}

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

    vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_quad_points];
    vector3 force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

    // Convert all area coordinate gauss points to cartesian
    for (int i = 0; i < num_tri_gauss_quad_points; i++) {
        f1->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
        f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
    }

    // Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
    // Also calculate energy whilst looping through face points
    scalar energy = 0.0;
    scalar mag_r, mag_ri,  mag_ri_2, mag_ri_4, mag_ri_6, mag_ri_7;
    scalar force_mag, vdw_fac_6;
    scalar vdw_r_eq_2 = vdw_r_eq * vdw_r_eq;
    scalar vdw_r_eq_4 = vdw_r_eq_2 * vdw_r_eq_2;
    scalar vdw_r_eq_6 = vdw_r_eq_4 * vdw_r_eq_2;
    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
        for(int l = k; l < num_tri_gauss_quad_points; l++) {
//          vector3 r =     {
//                                  minimum_image(p[k].x - q[l].x, box_size.x),
//                                  minimum_image(p[k].y - q[l].y, box_size.y),
//                                  minimum_image(p[k].z - q[l].z, box_size.z)
//                          };
            mag_r = sqrt( (p[k].x - q[l].x) * (p[k].x - q[l].x) + 
                          (p[k].y - q[l].y) * (p[k].y - q[l].y) +
                          (p[k].z - q[l].z) * (p[k].z - q[l].z) );
            mag_ri = 1./mag_r;
            mag_ri_2 = mag_ri * mag_ri;
            mag_ri_4 = mag_ri_2 * mag_ri_2;
            mag_ri_6 = mag_ri_4 * mag_ri_2;
            mag_ri_7 = mag_ri_6 * mag_ri;
            vdw_fac_6 = vdw_r_eq_6 * mag_ri_6;
            force_mag = 12 * mag_ri_7 * vdw_r_eq_6 * vdw_eps * (vdw_fac_6 - 1);
            energy += vdw_eps * vdw_fac_6 * (vdw_fac_6 - 2 );
            // force_mag *= -1;

            force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
            force_pair_matrix[k][l].y = force_mag * ((p[k].x - q[l].x) / mag_r);
            force_pair_matrix[k][l].z = force_mag * ((p[k].x - q[l].x) / mag_r);

            force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
            force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
            force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
        }
    }

    //			printf("YO\n");
    //			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
    //				for(int l = 0; l < num_tri_gauss_quad_points; l++) {
    //					printf("(%e %e %e) ", force_pair_matrix[k][l].x, force_pair_matrix[k][l].y, force_pair_matrix[k][l].z);
    //				}
    //				printf("\n");
    //			}

    scalar ApAq = f1->area * f2->area;
    energy *= ApAq;
    f1->add_bb_vdw_energy_to_record(energy, f2->daddy_blob->blob_index);
    f2->add_bb_vdw_energy_to_record(energy, f1->daddy_blob->blob_index);

    for (int j = 0; j < 3; j++) {
        vector3 force1 = {0, 0, 0}, force2 = {0, 0, 0};
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
        f1->add_bb_vdw_force_to_record(&force1, f2->daddy_blob->blob_index);
        //				printf("1:: %d %e %e %e\n", j, force1.x, force1.y, force1.z);

        force2.x *= ApAq;
        force2.y *= ApAq;
        force2.z *= ApAq;
        f2->add_force_to_node(j, &force2);
        f2->add_bb_vdw_force_to_record(&force1, f1->daddy_blob->blob_index);
        //				printf("2:: %d %e %e %e\n", j, force2.x, force2.y, force2.z);

        //				f2->add_force_to_node_atomic(j, &force);

        //				if(j == 0) {
        //					f2->add_force_to_node_atomic(0, &force);
        //				} else if(j == 1) {
        //					f2->add_force_to_node_atomic(1, &force);
        //				} else if(j == 2) {
        //					f2->add_force_to_node_atomic(2, &force);
        //				} else {
        //					printf("WTF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        //				}

        //				printf("Face 1: force on node %d = {%e %e %e}\n", j, force.x, force.y, force.z);
    }


    //			printf("\n\n\n\n");

    /*
                            for(int j = 0; j < 3; j++) {
                                    vector3 force = {0, 0, 0};
                                    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
                                            for(int l = 0; l < num_tri_gauss_quad_points; l++) {
                                                    scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];
                                                    force.x -= c * force_pair_matrix[k][l].x;
                                                    force.y -= c * force_pair_matrix[k][l].y;
                                                    force.z -= c * force_pair_matrix[k][l].z;
                                            }
                                    }
                                    force.x *= ApAq;
                                    force.y *= ApAq;
                                    force.z *= ApAq;
                                    f2->add_force_to_node(j, &force);
                                    printf("Face 2: force on node %d = {%e %e %e}\n", j, force.x, force.y, force.z);
                            }
                            printf("\n");

     */

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

    const int num_tri_gauss_quad_points = 3;

    const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] ={
        // Weight, eta1, eta2, eta3
        {0.333333333333333,
            {0.666666666666667, 0.166666666666667, 0.166666666666667}},
        {0.333333333333333,
            {0.166666666666667, 0.666666666666667, 0.166666666666667}},
        {0.333333333333333,
            {0.166666666666667, 0.166666666666667, 0.666666666666667}}

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
    f->add_xz_vdw_energy_to_record(energy);

    for (int j = 0; j < 3; j++) {
        vector3 force = {0, 0, 0};
        for (int k = 0; k < num_tri_gauss_quad_points; k++) {
            for (int l = 0; l < num_tri_gauss_quad_points; l++) {
                scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
                force.y += c * force_pair_matrix[k][l];
            }
        }
        force.y *= Asq;
        f->add_force_to_node(j, &force);
        f->add_xz_vdw_force_to_record(&force);
    }

}

scalar VdW_solver::distance2(vector3 *p, vector3 *q) {
    scalar dx = p->x - q->x, dy = p->y - q->y, dz = p->z - q->z;
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
