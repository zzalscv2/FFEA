#include "LJSteric_solver.h"

// This combined interaction is piecewise and as follows
//if r < 0, steric
//if 0 < r < rm, intermediate soft core potential
// if r > rm, lennard-jones

/** do_volumeExclusion calculates the force (and not the energy, yet) of two tetrahedra */
void LJSteric_solver::do_interaction(Face *f1, Face *f2){

    // First check two things (either of which results in not having to calculate anything):
    // Check that faces are facing each other, if not then they are not interacting
    if (dot(&f1->normal, &f2->normal) > ffea_const::zero) {
        return;
    }

    /* Robin suspects that this was leading to unstabilities...
     *  but this steric solver is more stable than the LJ one. */
    // Check that faces are in front of each other
    vector3 sep = {f2->centroid.x - f1->centroid.x, f2->centroid.y - f1->centroid.y, f2->centroid.z - f1->centroid.z};
    if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
        return;
    }

    //  Firstly, check that no nodes are shared,
    //     only in the case that faces belong to the same blob:
    if (f1->daddy_blob == f2->daddy_blob) {
      if (f1->n[3] == f2->n[3]) {
            return;
      }
      for (int i=0; i<4; i++) {
          int in_i = f1->n[i]->index;
          for (int j=0; j<4; j++) {
             if (f2->n[j]->index == in_i){
                 return;
             }
          }
      }
    }

    //  Then, check whether the tetrahedra intersect.
    if (f1->checkTetraIntersection(f2)) {
	    // If yes, do the vol-vol steric interactions
//cout << "Interacting!!" << endl;
//exit(0);
	    //   and get the direction of the force for f1:
	    /* TRIAL 2 */
	    grr3 force1, force2; //, n1_b;
	    vec3Vec3SubsToArr3(f2->n[3]->pos, f1->n[3]->pos, force2);
	    arr3Normalise<scalar,arr3>(force2); // that is the direction of the force for f1 (backwards).

	    /* TRIAL 1
	    arr3 force1, force2, n1_b;
	    vec3ResizeToArr3(ffea_const::mOne, f1->normal, n1_b);
	    vec3Arr3AddToArr3(f2->normal, n1_b, force1);
	    arr3Normalise<scalar,arr3>(force1); // that is the direction of the force for f1 (backwards).
	    */

	    //////////////////////////////////////////////
	    // One more check ////// One more check /////
	    /*
	    arr3 inwards;
	    vec3Vec3SubsToArr3(f1->n[3]->pos, f1->centroid, inwards);
	    if (arr3arr3DotProduct<scalar,arr3>(inwards, force1) < ffea_const::zero) {
	      return;
	    }
	    */
	    // end checking //////////////////////////////

	    // scalar vol_f = 1
	    // Finally, get the intersection volume:
	    // scalar vol = f1->getTetraIntersectionVolume(f2);
	    geoscalar vol, dVdr;
	    f1->getTetraIntersectionVolumeAndGradient(f2,force2,vol,dVdr);

	    vol *= steric_factor;
	    dVdr *= steric_factor;

	    // Energy is proportional to the volume of interaction:
	    //f1->add_bb_vdw_energy_to_record(vol, f2->daddy_blob->blob_index);
	    //f2->add_bb_vdw_energy_to_record(vol, f1->daddy_blob->blob_index);

	    // Store the measurement
	    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;

	    // arr3Resize(vol, force1);  // the provious volume force
	    // Force is proportional to the gradient of this volume:
	    arr3Resize<scalar,arr3>(scalar(-dVdr), force2);
	    arr3Resize2<scalar,arr3>(ffea_const::mOne, force2, force1);

	    for (int j = 0; j < 3; j++) {
	      f1->add_force_to_node(j, force1);
	      f1->add_bb_vdw_force_to_record(force1, f2->daddy_blob->blob_index);
	      f2->add_force_to_node(j, force2);
	      f2->add_bb_vdw_force_to_record(force2, f1->daddy_blob->blob_index);
	    }

    } else {

      // If not, do lennard-jones type interactions

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
	    scalar force_mag, vdw_fac, vdw_fac_2, vdw_fac_6;

	    scalar vdw_r_eq_2 = vdw_r_eq * vdw_r_eq;
	    scalar vdw_r_eq_4 = vdw_r_eq_2 * vdw_r_eq_2;
	    scalar vdw_r_eq_6 = vdw_r_eq_4 * vdw_r_eq_2;

	    scalar vdw_r_eqi = 1.0 / vdw_r_eq;
	    scalar vdw_r_eqi2 = vdw_r_eqi * vdw_r_eqi;
	    // scalar vdw_r_eqi3 = vdw_r_eqi * vdw_r_eqi2;

	    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
		for(int l = k; l < num_tri_gauss_quad_points; l++) {

		    mag_r = sqrt( (p[k].x - q[l].x) * (p[k].x - q[l].x) +
		                  (p[k].y - q[l].y) * (p[k].y - q[l].y) +
		                  (p[k].z - q[l].z) * (p[k].z - q[l].z) );

		    // Both parts need this
		    mag_ri = 1./mag_r;

		    if(mag_r < vdw_r_eq) {

			// Intermediatey stuff

			vdw_fac = mag_r * vdw_r_eqi;
			vdw_fac_2 = vdw_fac * vdw_fac;
			force_mag = 6 * vdw_eps * mag_ri * vdw_fac_2 * (vdw_fac - 1);
			energy += vdw_eps * vdw_fac_2 * (2 * vdw_fac - 3);

			force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
			force_pair_matrix[k][l].y = force_mag * ((p[k].x - q[l].x) / mag_r);
			force_pair_matrix[k][l].z = force_mag * ((p[k].x - q[l].x) / mag_r);

			force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
			force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
			force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;

		    } else {
			    mag_ri_2 = mag_ri * mag_ri;
			    mag_ri_4 = mag_ri_2 * mag_ri_2;
			    mag_ri_6 = mag_ri_4 * mag_ri_2;
			    mag_ri_7 = mag_ri_6 * mag_ri;
			    vdw_fac_6 = vdw_r_eq_6 * mag_ri_6;
			    force_mag = 12 * mag_ri_7 * vdw_r_eq_6 * vdw_eps * (vdw_fac_6 - 1);
			    energy += vdw_eps * vdw_fac_6 * (vdw_fac_6 - 2 );

			    force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
			    force_pair_matrix[k][l].y = force_mag * ((p[k].x - q[l].x) / mag_r);
			    force_pair_matrix[k][l].z = force_mag * ((p[k].x - q[l].x) / mag_r);

			    force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
			    force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
			    force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
		    }
		}
	    }


	    scalar ApAq = f1->area * f2->area;
	    energy *= ApAq;
	    //f1->add_bb_vdw_energy_to_record(energy, f2->daddy_blob->blob_index);
	    //f2->add_bb_vdw_energy_to_record(energy, f1->daddy_blob->blob_index);

	    // Store the measurement
	    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += energy;

	    for (int j = 0; j < 3; j++) {
		vector3 force1 = {0, 0, 0}, force2 = {0, 0, 0};
		for (int k = 0; k < num_tri_gauss_quad_points; k++) {
		    for (int l = 0; l < num_tri_gauss_quad_points; l++) {
		        scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
		        scalar d = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];

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

		force2.x *= ApAq;
		force2.y *= ApAq;
		force2.z *= ApAq;
		f2->add_force_to_node(j, &force2);
		f2->add_bb_vdw_force_to_record(&force2, f1->daddy_blob->blob_index);
	    }

    }

    return;

}
    /**Calculates LJSteric forces modified with periodic boundary correction in distance calculation*/
void LJSteric_solver::do_interaction(Face *f1, Face *f2, scalar *blob_corr){


    int f1_daddy_blob_index = f1->daddy_blob->blob_index;
    int f2_daddy_blob_index = f2->daddy_blob->blob_index;
    // First check two things (either of which results in not having to calculate anything):
    // Check that faces are facing each other, if not then they are not interacting


    if (dot(&f1->normal, &f2->normal) > ffea_const::zero) {
        return;
    }

    /* Robin suspects that this was leading to unstabilities...
     *  but this steric solver is more stable than the LJ one. */
    // Check that faces are in front of each other (with periodicity modification)
    vector3 sep = {f2->centroid.x - f1->centroid.x-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3],f2->centroid.y - f1->centroid.y-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3+1],f2->centroid.z - f1->centroid.z-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3+2]};
    if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
        return;
    }

    //  Firstly, check that no nodes are shared,
    //     only in the case that faces belong to the same blob:
    if (f1->daddy_blob == f2->daddy_blob) {
      if (f1->n[3] == f2->n[3]) {
            return;
      }
      for (int i=0; i<4; i++) {
          int in_i = f1->n[i]->index;
          for (int j=0; j<4; j++) {
             if (f2->n[j]->index == in_i){
                 return;
             }
          }
      }
    }

    //  Then, check whether the tetrahedra intersect.(with periodicity modification)
    if (f1->checkTetraIntersection(f2, blob_corr, f1_daddy_blob_index, f2_daddy_blob_index)) {
	    // If yes, do the vol-vol steric interactions
//cout << "Interacting!!" << endl;
//exit(0);
	    //   and get the direction of the force for f1 (with periodicity modification):
	    /* TRIAL 2 */
	    grr3 force1, force2; //, n1_b;
	    f1->vec3Vec3SubsToArr3Mod(f2, force2, blob_corr, f1_daddy_blob_index, f2_daddy_blob_index);
	    arr3Normalise<scalar,arr3>(force2); // that is the direction of the force for f1 (backwards).

	    /* TRIAL 1
	    arr3 force1, force2, n1_b;
	    vec3ResizeToArr3(ffea_const::mOne, f1->normal, n1_b);
	    vec3Arr3AddToArr3(f2->normal, n1_b, force1);
	    arr3Normalise<scalar,arr3>(force1); // that is the direction of the force for f1 (backwards).
	    */

	    //////////////////////////////////////////////
	    // One more check ////// One more check /////
	    /*
	    arr3 inwards;
	    vec3Vec3SubsToArr3(f1->n[3]->pos, f1->centroid, inwards);
	    if (arr3arr3DotProduct<scalar,arr3>(inwards, force1) < ffea_const::zero) {
	      return;
	    }
	    */
	    // end checking //////////////////////////////

	    // scalar vol_f = 1
	    // Finally, get the intersection volume (with periodicity modification):
	    // scalar vol = f1->getTetraIntersectionVolume(f2);
	    geoscalar vol, dVdr;
	    f1->getTetraIntersectionVolumeAndGradient(f2,force2,vol,dVdr,blob_corr, f1_daddy_blob_index,f2_daddy_blob_index);

	    vol *= steric_factor;
	    dVdr *= steric_factor;

	    // Energy is proportional to the volume of interaction:
	    //f1->add_bb_vdw_energy_to_record(vol, f2->daddy_blob->blob_index);
	    //f2->add_bb_vdw_energy_to_record(vol, f1->daddy_blob->blob_index);

	    // Store the measurement
	    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;

	    // arr3Resize(vol, force1);  // the provious volume force
	    // Force is proportional to the gradient of this volume:
	    arr3Resize<scalar,arr3>(scalar(-dVdr), force2);
	    arr3Resize2<scalar,arr3>(ffea_const::mOne, force2, force1);

	    for (int j = 0; j < 3; j++) {
	      f1->add_force_to_node(j, force1);
	      f1->add_bb_vdw_force_to_record(force1, f2->daddy_blob->blob_index);
	      f2->add_force_to_node(j, force2);
	      f2->add_bb_vdw_force_to_record(force2, f1->daddy_blob->blob_index);
	    }

    } else {

      // If not, do lennard-jones type interactions

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

	    };

	    vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_quad_points];
	    vector3 force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

	    // Convert all area coordinate gauss points to cartesian (with periodicity modification on one point)
	    for (int i = 0; i < num_tri_gauss_quad_points; i++) {
		f1->barycentric_calc_point_f2(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i],blob_corr,f1_daddy_blob_index,f2_daddy_blob_index);
		f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
	    }

	    // Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
	    // Also calculate energy whilst looping through face points
	    scalar energy = 0.0;
	    scalar mag_r, mag_ri,  mag_ri_2, mag_ri_4, mag_ri_6, mag_ri_7;
	    scalar force_mag, vdw_fac, vdw_fac_2, vdw_fac_6;

	    scalar vdw_r_eq_2 = vdw_r_eq * vdw_r_eq;
	    scalar vdw_r_eq_4 = vdw_r_eq_2 * vdw_r_eq_2;
	    scalar vdw_r_eq_6 = vdw_r_eq_4 * vdw_r_eq_2;

	    scalar vdw_r_eqi = 1.0 / vdw_r_eq;
	    scalar vdw_r_eqi2 = vdw_r_eqi * vdw_r_eqi;
	    // scalar vdw_r_eqi3 = vdw_r_eqi * vdw_r_eqi2;

	    for(int k = 0; k < num_tri_gauss_quad_points; k++) {
		for(int l = k; l < num_tri_gauss_quad_points; l++) {

		    mag_r = sqrt( (p[k].x - q[l].x) * (p[k].x - q[l].x) +
		                  (p[k].y - q[l].y) * (p[k].y - q[l].y) +
		                  (p[k].z - q[l].z) * (p[k].z - q[l].z) );

		    // Both parts need this
		    mag_ri = 1./mag_r;

		    if(mag_r < vdw_r_eq) {

			// Intermediatey stuff

			vdw_fac = mag_r * vdw_r_eqi;
			vdw_fac_2 = vdw_fac * vdw_fac;
			force_mag = 6 * vdw_eps * mag_ri * vdw_fac_2 * (vdw_fac - 1);
			energy += vdw_eps * vdw_fac_2 * (2 * vdw_fac - 3);

			force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
			force_pair_matrix[k][l].y = force_mag * ((p[k].x - q[l].x) / mag_r);
			force_pair_matrix[k][l].z = force_mag * ((p[k].x - q[l].x) / mag_r);

			force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
			force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
			force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;

		    } else {
			    mag_ri_2 = mag_ri * mag_ri;
			    mag_ri_4 = mag_ri_2 * mag_ri_2;
			    mag_ri_6 = mag_ri_4 * mag_ri_2;
			    mag_ri_7 = mag_ri_6 * mag_ri;
			    vdw_fac_6 = vdw_r_eq_6 * mag_ri_6;
			    force_mag = 12 * mag_ri_7 * vdw_r_eq_6 * vdw_eps * (vdw_fac_6 - 1);
			    energy += vdw_eps * vdw_fac_6 * (vdw_fac_6 - 2 );

			    force_pair_matrix[k][l].x = force_mag * ((p[k].x - q[l].x) / mag_r);
			    force_pair_matrix[k][l].y = force_mag * ((p[k].x - q[l].x) / mag_r);
			    force_pair_matrix[k][l].z = force_mag * ((p[k].x - q[l].x) / mag_r);

			    force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
			    force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
			    force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
		    }
		}
	    }


	    scalar ApAq = f1->area * f2->area;
	    energy *= ApAq;
	    //f1->add_bb_vdw_energy_to_record(energy, f2->daddy_blob->blob_index);
	    //f2->add_bb_vdw_energy_to_record(energy, f1->daddy_blob->blob_index);

	    // Store the measurement
	    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += energy;

	    for (int j = 0; j < 3; j++) {
		vector3 force1 = {0, 0, 0}, force2 = {0, 0, 0};
		for (int k = 0; k < num_tri_gauss_quad_points; k++) {
		    for (int l = 0; l < num_tri_gauss_quad_points; l++) {
		        scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
		        scalar d = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];

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

		force2.x *= ApAq;
		force2.y *= ApAq;
		force2.z *= ApAq;
		f2->add_force_to_node(j, &force2);
		f2->add_bb_vdw_force_to_record(&force2, f1->daddy_blob->blob_index);
	    }

    }

    return;

}

