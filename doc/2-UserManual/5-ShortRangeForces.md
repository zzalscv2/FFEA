

Short Range Forces {#shortRange}
================================



Overview {#srOverview}
======================


Two different short range forces have been implemented in FFEA. They are 
 mutually exclusive, so you need to choose which one to use. In any case, 
 you need to set:

     < calc_vdw = 1 > 

to turn them on. In addition, these short ranged forces are acting between 
 active faces. Active faces will be read from the file:

     < vdw = my-system.vdw >

where "my-system" is whatever the file you have, and this field is 
 included into the ` <conformation> ` block in the input .ffea file.
 The ` vdw ` file is a text file that starts with the lines:

     ffea vdw file
     num_faces Number-Of-Active-Faces
     vdw params:

where "Number-Of-Active-Faces" is the number of faces that your system has 
 pointing outwards, and every line can have values ranging from -1, 
 corresponding to inactive, up to 7 for the rest of faces types. 
 This file should have been automatically generated when
  [configuring the system](\ref makeffeablob), 
 with all the faces set up to "inactive", 
 and it can be easily configured using the [FFEA viewer](\ref FFEAviewer).

The next thing you need is to put a box in your system. This means giving values for
 ` es_N_x `, ` es_N_y ` and ` es_N_z `, as well as for ` kappa `. See the
 [keyword reference](ref paramBlock) to get information on the values that these
 fields should have. 

where es_N_xyz fields are how many voxels your simulation box has in each direction. The size of each voxel is es_h * 1/kappa, where kappa is the inverse debye length. So basically, each voxel is of size how many debye lengths should things communicate over. I usually use three. Your simulation box should contain your entire molecule set of course which is where you decide es_N_xyz, as in how many voxels do you need in each direction. The box defaults to periodic boundary conditions.

    




Lennard-Jones potential {#ljPotential}
======================================

The well established 6-12 Lennard-Jones potential is used if you 




Steric potential {#sPotential}
==============================



