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

    // //  Working version for F = k*dV/dr // //
    if (! f1->checkTetraIntersection(f2)) return;
    geoscalar vol;
    grr3 dVdr;
    grr4 phi1, phi2;

    if (!f1->getTetraIntersectionVolumeTotalGradientAndShapeFunctions(f2, dVdr, vol, phi1, phi2)) return;

    vol *= steric_factor;

    // Force is proportional to the gradient, i. e.:
    arr3Resize<geoscalar,grr3>(steric_factor, dVdr);

    grr3 ftmp1, ftmp2;
    #pragma omp critical 
    { 
    // Store the measurement
    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;
    // Finally, apply the force onto the nodes:
    for (int j = 0; j < 4; j++) {
      arr3Resize2<geoscalar,grr3>(phi1[j], dVdr, ftmp1);
      f1->add_force_to_node(j, ftmp1);
      f1->add_bb_vdw_force_to_record(ftmp1, f2->daddy_blob->blob_index);

      arr3Resize2<geoscalar,grr3>(ffea_const::mOne*phi2[j], dVdr, ftmp2);
      f2->add_force_to_node(j, ftmp2);
      f2->add_bb_vdw_force_to_record(ftmp2, f1->daddy_blob->blob_index);
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


}

    /**Calculates Steric forces modified with periodic boundary correction in distance calculation*/
void Steric_solver::do_interaction(Face *f1, Face *f2, scalar * blob_corr){


    //printf("Correct interaction triggered\n");
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

    if (! f1->checkTetraIntersection(f2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index)) return;
    geoscalar vol;
    grr3 dVdr;
    grr4 phi1, phi2;

    if (!f1->getTetraIntersectionVolumeTotalGradientAndShapeFunctions(f2, dVdr, vol, phi1, phi2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index)) return;

    vol *= steric_factor;

    // Force is proportional to the gradient, i. e.:
    arr3Resize<geoscalar,grr3>(steric_factor, dVdr);

    grr3 ftmp1, ftmp2;
    // Store the measurement
    #pragma omp critical 
    { 
    fieldenergy[f1->daddy_blob->blob_index][f2->daddy_blob->blob_index] += vol;
    // Finally, apply the force onto the nodes:
    for (int j = 0; j < 4; j++) {
      arr3Resize2<geoscalar,grr3>(phi1[j], dVdr, ftmp1);
      f1->add_force_to_node(j, ftmp1);
      f1->add_bb_vdw_force_to_record(ftmp1, f2->daddy_blob->blob_index);

      arr3Resize2<geoscalar,grr3>(ffea_const::mOne*phi2[j], dVdr, ftmp2);
      f2->add_force_to_node(j, ftmp2);
      f2->add_bb_vdw_force_to_record(ftmp2, f1->daddy_blob->blob_index);
    }
    }






    /* //  Then, check whether the tetrahedra intersect(with periodicity modification),
    //    and if so, get the volume:
    scalar vol = f1->checkTetraIntersectionAndGetVolume(f2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index);
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
      cm2[0] += f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3];
      cm2[1] += f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3+1];
      cm2[2] += f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3+2];
    }
    arr3Resize<scalar,arr3>(0.25,cm1);
    arr3Resize<scalar,arr3>(0.25,cm2);
    arr3arr3Substract<scalar,arr3>(cm2, cm1, force2);
    //printf("**********\n Blob %d to Blob %d\n face %d to face %d\n corr.x is %f\ndist in x is %f\ndist in y is %f\ndist in z is %f\n",f1_daddy_blob_index,f2_daddy_blob_index,f1->index, f2->index,blob_corr[f1_daddy_blob_index*(num_blobs)*3 + f2_daddy_blob_index*3],force2[0],force2[1],force2[2]);
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
    }*/

}
