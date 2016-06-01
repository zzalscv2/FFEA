
Input and Output files {#ioFiles}
======================

Input files {#iFiles}
=====================

Main runner input file: .ffea  
------------------------------
 Corresponding to the main input file for the FFEA runner, this file is widely described 
  [here](\ref ffea_iFile), and [here](\ref keywordReference)



Nodes file: .node
-----------------

Contains a list of whitespace separated coordinates specifiying the position
 of each node in the blob. The file starts with a header of four lines, and 
 is followed by two separated lists with the interior and exterior nodes:


    ffea node file
    num_nodes 305
    num_surface_nodes 242
    num_interior_nodes 63
    surface nodes:
    x0 y0 z0
    x1 y1 z1
    ...
    x241 y241 z241
    interior nodes:
    x242 y242 z242
    x243 y243 z243
    ...
    x304 y304 z304



Topology file: .top
-------------------

Contains a list a whitespace separated indices specifiying the node indices
contained within each element. There are ten indices per element as they are second
order elements (every tetrahedron has a node at each corner 
 and at the midpoint of each edge). Separated into internal (no surface faces) 
 and external (surface faces) elements for use in surface-surface interactions.


    ffea topology file
    num_elements 130
    num_surface elements 120
    num_interior elements 10
    surface elements:
    n00 n01 n02 n03 n04 n05 n06 n07 n08 n09
    n10 n11 n12 n13 n14 n15 n16 n17 n18 n19
    ...
    n110 n111 n112 n113 n114 n115 n116 n117 n118 n119
    interior elements:
    n120 n121 n122 n123 n124 n125 n126 n127 n128 n129



Surface file: .surf
-------------------

Contains a list of whitespace separated indices specifying the indices of the
nodes on each surface face and the index of the containing element.


    ffea surface file
    num_surface_faces 480
    faces:
    el0 n00 n01 n02
    el1 n10 n11 n12
    ... 
    el480 n4800 n4801 n4802



Material file: .mat
-------------------

Contains a list of whitespace separated material parameters for each element.


    ffea material params file
    num elements 120
    density0 shear_visc0 bulk_visc0 shear_modulus0 bulk_modulus0 dielectric_const0
    density1 shear_visc1 bulk_visc1 shear_modulus1 bulk_modulus1 dielectric_const1
    ...
    density119 shear_visc119 bulk_visc119 shear_modulus119 bulk_modulus119 dielectric_const119
    


Stokes radii file: .stokes
--------------------------

Contains a list of effective hydrodynamic radii for each node in the system
for the calculation of external drag.

    ffea stokes radii file
    num_nodes 305
    radius0
    radius1
    ... 
    radius 304


van der Waals file: .vdw
-------------------------
Contains a list of integers, ranging from -1 to 6, describing the type of vdw interaction
 this face will undergo, as described by the Lennard-Jones files, where:
   * -1 - stands for inactive faces.
   * 0 to 6 - specifies the type of an active face, interacting 
              with parameters defined in .lj file.


    ffea vdw file
    num_faces 480
    vdw_params:
    type0
    type1
    ...
    type479


Lennard-Jones file: .lj
-----------------------

Contains a matrix of parameter pairs for the Lennard-Jones interactions 
 between the different face types specified in the .vdw file. 
 A description of this file can be found in the [Lennard-Jones](\ref ljPotential) section. 

<!-- In the following example:


    ffea vdw forcefield params file
    num vdw face types 7
    (1e+15, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)


face type0 interacts with face type0 through a Lennard-Jones potential
 that has a well depth of \f$\epsilon_{0,0} = 10^{15} J/m^2\f$ and an 
 equilibrium distance \f$\sigma_{0,0} = 10^{-9} m\f$, 
 whereas face type 0 interacts with face type 1 
 (and in fact, all other interaction types) with a Lennard-Jones potential
 with parameters \f$\epsilon_{0,1} = 10^{12} J/m^2\f$ and 
  \f$\sigma_{0,1} = 10^{-9} m\f$. Therefore, this matrix should be symmetric. 
--> 


Pre-computed potentials: .pot and .force  {#potfile}
----------------------------------------
These files contain values for the potential and force in a two column format:


    0 1606.41515966667
    0.002 1595.931897
    0.004 1580.257003
    0.006 1564.662109


where the first column give positions and the second column has energies and forces
 respectively. The runner input file will expect a conversion factor to \f$m\f$ and 
 \f$J\f$ (and \f$Jm\f$) for them (see the [pre-computed potentials](\ref fm_inputfile)
 for more details). 


Springs file: .springs
----------------------

This file contains a set of springs joining pairs of nodes. After an introductory header, 
 every line specifies the details of a single spring through:
 equilibrium distance (in meters), potential well depth (in ` N/m `), 
 blob i, blob j, conformation k, conformation l, node m, node n. 
Specifically, a springs file with 3 springs could look like:


    ffea springs file
    num_springs 3
    springs:
    1e-2 5e-9 0 1 0 0 40 40
    1e-2 5e-9 0 1 0 0 57 57
    1e-2 3e-9 0 1 0 0 66 83





States file
-----------


Rates file
----------


Map file: .map
--------------




Output files {#oFiles}
======================

Trajectory file: .trj
---------------------
Contains a list of `*` separated frames specifying the structure of the blob at
each outputted timestep. The frames contain a list of blobs which themselves
contain a list of node positions, velocities and forces. If the motion state is
STATIC, blob doesn’t move, so frame is not outputted to save space, and improve 
the performance. Only the positions are used in the viewer.
 In between each frame is a small segment showing whether any conformational 
 changes have occured. As such,
each frame can be a different length, which is why the blob and conformation
sizes are specified at the beginning of the file.


    FFEA trajectory file

    Initialisation:
    Number of Blobs 2
    Number of Conformations 1 1
    Blob 0: Conformation 0 Nodes 305
    Blob 1: Conformation 0 Nodes 305

    *
    Blob 0, Conformation 0, step 0
    DYNAMIC
    x0 y0 z0 velx0 vely0 velz0 phi0 forcex0 forcey0 forcez0
    x1 y1 z1 velx1 vely1 velz1 phi1 forcex1 forcey1 forcez1
    ... 
    Blob 1, Conformation 0, step 0 DYNAMIC x0 y0 z0 velx0 vely0 velz0 phi0
    forcex0 forcey0 forcez0
    x1 y1 z1 velx1 vely1 velz1 phi1 forcex1 forcey1 forcez1
    ...
    Conformation Changes:
    Blob 0: Conformation 0 -> Conformation 1
    Blob 1: Conformation 0 -> Conformation 0
    *
    Blob 0, Conformation 1, step 1001
    ... 




Measurement files: .out
------------------------
 
Contains a list of the relevant system properties of each outputted timestep.
Two seperate types of file for internal and external measurements.


### Blob measurements: .out ### 

There is a file of this kind for every blob, containing a 
 list of internal properties for each outputted timestep. Specifically,
 all energies, centres of mass, momenta, rmsd and vdw interactions with the
 simulation box (if sticky wall xz = 1) are recorded in these files. The format
 is:


    # step — KE — PE — CoM x — CoM y — CoM z — L x — L y — L z — rmsd — vdw area 0 surface — vdw force 0 surface — vdw energy 0 surface
    ... 


### World measurements: .out ### 

This file records the Lennard-jones 
  interaction values between each pair of blobs in the whole system along the trajectory:


    # step — vdw area 0 1 — vdw force 0 1 — vdw energy 0 1
.
.

