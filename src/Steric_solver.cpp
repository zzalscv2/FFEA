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

    if (! f1->checkTetraIntersection(f2)) return;
    do_steric_interaction(f1, f2); 
    return;

}

    /**Calculates Steric forces modified with periodic boundary correction in distance calculation*/
void Steric_solver::do_interaction(Face *f1, Face *f2, scalar * blob_corr){

    int f1_daddy_blob_index = f1->daddy_blob->blob_index;
    int f2_daddy_blob_index = f2->daddy_blob->blob_index;

    if (! f1->checkTetraIntersection(f2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index)) return;
    geoscalar vol;
    grr3 dVdr;
    grr4 phi1, phi2;

    if (!f1->getTetraIntersectionVolumeTotalGradientAndShapeFunctions(f2, steric_dr, dVdr, vol, phi1, phi2,blob_corr,f1_daddy_blob_index, f2_daddy_blob_index)) return;

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
      // f1->add_bb_vdw_force_to_record(ftmp1, f2->daddy_blob->blob_index); // DEPRECATED

      arr3Resize2<geoscalar,grr3>(ffea_const::mOne*phi2[j], dVdr, ftmp2);
      f2->add_force_to_node(j, ftmp2);
      // f2->add_bb_vdw_force_to_record(ftmp2, f1->daddy_blob->blob_index); // DEPRECATED
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
