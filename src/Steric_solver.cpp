#include "Steric_solver.h"


/** do_volumeExclusion calculates the force (and not the energy, yet) of two tetrahedra */
void Steric_solver::do_interaction(Face *f1, Face *f2){ 

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
    if (!f1->checkTetraIntersection(f2)) {
	return;
 /*   #pragma omp critical
    {
	fprintf(stderr, "%d %d Nointersect", f1->index, f2->index);
	for(int i = 0; i < 4; ++i) {
		fprintf(stderr, " %6.3f %6.3f %6.3f", f1->n[i]->pos.x, f1->n[i]->pos.y, f1->n[i]->pos.z);
		fprintf(stderr, "  ");
	}
	for(int i = 0; i < 4; ++i) {
		fprintf(stderr, " %6.3f %6.3f %6.3f", f2->n[i]->pos.x, f2->n[i]->pos.y, f2->n[i]->pos.z);
		fprintf(stderr, "  ");	
	}
	fprintf(stderr, "\n");
    }*/
    } else {
  /*  #pragma omp critical
    {
	fprintf(stderr, "%d %d Intersect", f1->index, f2->index);
	for(int i = 0; i < 4; ++i) {
		fprintf(stderr, " %6.3f %6.3f %6.3f", f1->n[i]->pos.x, f1->n[i]->pos.y, f1->n[i]->pos.z);
		fprintf(stderr, "  ");
	}
	for(int i = 0; i < 4; ++i) {
		fprintf(stderr, " %6.3f %6.3f %6.3f", f2->n[i]->pos.x, f2->n[i]->pos.y, f2->n[i]->pos.z);
		fprintf(stderr, "  ");	
	}
	fprintf(stderr, "\n");
    }*/
    }

   /* fprintf(stderr, "Interacting Faces:\n");
    fprintf(stderr, "%d %d %d %d\n", f1->index, f1->n[0]->index, f1->n[1]->index, f1->n[2]->index);
    fprintf(stderr, "%d %d %d %d\n", f2->index, f2->n[0]->index, f2->n[1]->index, f2->n[2]->index);
    f1->print_centroid();
    f1->print_nodes();
    f2->print_centroid();
    f2->print_nodes();*/
  //  exit(0);
    //   and get the direction of the force for f1:
    /* TRIAL 2 */
    arr3 force1, force2; //, n1_b; 
    vec3Vec3SubsToArr3(f1->n[3]->pos, f2->n[3]->pos, force1);
    arr3Normalise<scalar,arr3>(force1); // that is the direction of the force for f1 (backwards). 

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
    geoscalar vol, area; 
    f1->getTetraIntersectionVolumeAndArea(f2,vol,area);
    //fprintf(stderr, "Vol = %f\nArea = %f\n", vol, area);
    //exit(0);
    area *= steric_factor;
    vol *= steric_factor; 

    // Energy is proportional to the volume of interaction: 
//    f1->add_bb_vdw_energy_to_record(vol, f2->daddy_blob->blob_index);
//    f2->add_bb_vdw_energy_to_record(vol, f1->daddy_blob->blob_index);

    // Store the measurement
    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;

    // arr3Resize(vol, force1);  // the provious volume force 
    // Force is proportional to the surface area of this volume: 
    arr3Resize(scalar(area), force1); 
    arr3Resize2(ffea_const::mOne, force1, force2);
    
    for (int j = 0; j < 3; j++) {
      f1->add_force_to_node(j, force1); 
      f1->add_bb_vdw_force_to_record(force1, f2->daddy_blob->blob_index);
      f2->add_force_to_node(j, force2); 
      f2->add_bb_vdw_force_to_record(force2, f1->daddy_blob->blob_index);
    } 


}

