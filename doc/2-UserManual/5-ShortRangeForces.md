

Short Range Forces {#shortRange}
================================



Overview {#srOverview}
======================


Two different short range forces have been implemented in FFEA. They are 
 mutually exclusive, so you need to choose which one to use. In any case, 
 you need to set:

     <calc_vdw = 1> 

to enable it. In addition, these short ranged forces are acting between 
 active faces. Active faces will be read from the file:

     <vdw = my-system.vdw >

where "my-system" is whatever the file you have, and this field is 
 included into the ` <conformation> ` block. That file is a text file 
 that starts with the lines:

     ffea vdw file
     num_faces Number-Of-Active-Faces
     vdw params:

where "Number-Of-Active-Faces" is the number of faces that your system has 
 pointing outwards, and every line can have values ranging from -1, 
 corresponding to inactive, to 7 for the rest of faces types. 
 This file should have been automatically generated when
  [configuring the system](\ref makeffeablob), 
 with all the faces set up to "inactive", 
 and it can be easily configured using the [FFEA viewer](\ref FFEAviewer).





Lennard-Jones potential {#ljPotential}
======================================

The well established 6-12 Lennard-Jones potential is used if you 




Steric potential {#sPotential}
==============================



