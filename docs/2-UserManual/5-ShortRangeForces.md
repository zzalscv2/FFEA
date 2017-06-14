

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
 included in the `<blob>` block (or in the ` <conformation> ` block, 
 if multiple conformations are defined) in the input .ffea file.
 The ` vdw ` file is a text file that starts with the lines:

     ffea vdw file
     num_faces Number-Of-Active-Faces
     vdw params:

where "Number-Of-Active-Faces" is the number of faces that your system has 
 pointing outwards, and every line can have values ranging from -1, 
 corresponding to **inactive**, up to 7 for the rest of faces types. 
 This file should have been automatically generated when
  [configuring the system](\ref voltoffeatut), 
 with all the faces set up to "inactive". There are a number of ways of
 configuring this file with different values for different files, either 
 geometrically or through direct visualisation.
 The most obvious approach may be to use the [FFEA viewer](\ref FFEAviewertut)
 showing the ` Face Indices ` in ` Show Indices `. However, reading the 
 interesting nodes from the PyMOL viewer window may be quite tricky. Thus,
 one may prefer 
 to ` Add Atoms ` ` Onto Nodes `, and then select the nodes of interest by
 clicking on the screen. The selection may be printed typing:

     print cmd.get_pdbstr("sele")

 being "sele" the selected atoms. The residue number printed out will 
 correspond to the internal face index (starting with 0), which corresponds
 to the order of the faces in the ` vdw ` file. 


Another keyword ` inc_self_vdw ` will determine whether interactions due 
 to interacting faces within the same blob are taken into account (default, 1) or not (0). 

The next thing you need is to put a box in your system. This means giving values for
 ` es_h `, ` es_N_x `, ` es_N_y ` and ` es_N_z `, as well as for ` kappa `. 
 Even if [keyword reference](\ref paramBlock) define these keywords, we introduce 
 them here, because FFEA will only compute face-face interactions that fall into 
 the same or adjacent cells. The number of voxels that the simulation box has in each 
 direction is defined through ` es_N_x `, ` es_N_y ` and ` es_N_z `, and the size of each voxel 
 is ` es_h ` / ` kappa `, where ` kappa ` is the inverse Debye length.

    




Lennard-Jones potential {#ljPotential}
======================================

In the case of setting:

    < vdw_type = lennard-jones > 

the well known 6-12 Lennard-Jones potential:
\f[

  U_{i,j}(r) = \epsilon_{i,j} \left[ \left( \frac{\sigma_{i,j}}{r} \right)^{12} - 2 \left( \frac{\sigma_{i,j}}{r} \right)^6 \right]

\f]

is used to measure surface-surface interactions between all the possible set of 
 face pairs, where \f$\sigma_{i,j}\f$ is the equilibrium separation distance, and 
 \f$\epsilon_{i,j}\f$ is the depth of the energy well at \f$r = r_{i,j}\f$ 
 per surface unit for two faces interacting of faces `i` and `j`. 
 A `lj` file containing the different \f$\sigma_{ij}\f$ and \f$\epsilon_{i,j}\f$ 
 interacting parameters for all the possible face pairs is needed as input. 
 This file, has to be given in the .ffea as:
   
     < vdw_forcefield_params = your-file.lj >

and has the format:


     ffea vdw forcefield params file
     num_vdw_face_types 7
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+15, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+15, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)
     (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10) (1e+13, 5e-10)

where the matrix provides the values (\f$\epsilon_{i,j}\f$, \f$\sigma_{i,j}\f$) in 
  \f$J/m^2\f$ and \f$m\f$, respectively. Because of being a surface-surface interaction, 
 the larger the faces interacting the stronger the interaction. 





Steric potential {#sPotential}
==============================
 In the case of setting:

      < vdw_type = steric >

 a steric repulsion that is proportional to the overlapping of the tetrahedra 
  of the interacting faces will be calculated. More specifically, the repulsive 
  energy will be proportional to the volume of the overlapping tetrahedra, and 
  the repulsive force proportional to the gradient of this volume in the direction
  of the line joining the node in every tetrahedra that is not part of the interacting face.
 The magnitude of this interaction can be modulated through a single parameter, 
  
      < vdw_steric_factor = F >

 where ` F ` is a value to be provided. In the case of being negative, the 
  user will receive a warning. 

Combination potential {#cPotential}
==============================
 In the case of setting:

      < vdw_type = ljsteric >

 a piecewise combination of the steric potential and the lennard-jones potentials is used.
  Because hard-core surface - surface Lennard-Jones repulsion is much less stable than the softer volume - volume
  steric repulsion, we use the following protocol:

\f[  U_{i,j}(r) = \left\{
\begin{array}{ll}
      steric & r\leq 0 \\
      \epsilon_{i,j} \left[ 2\left(\frac{r}{\sigma_{i,j}}\right)^{3} - 3\left(\frac{r}{\sigma_{i,j}} \right)^2 \right]  & 0 < r\leq \sigma_{i,j} \\
      \epsilon_{i,j} \left[ \left(\frac{\sigma_{i,j}}{r}\right)^{12} - 2\left(\frac{\sigma_{i,j}}{r} \right)^6 \right]& r > \sigma_{i,j} \\
\end{array} 
\right. \f]

 Very short range interactions are dealt with using the volume-volume steric interactions, and `long range' using standard 
  Lennard-Jones interactions. The intermediate region uses an interpolated function, for which the function itself and it's first derivative
  are both continuous at the boundaries, and for \f$r = \sigma_{i,j}\f$, we still find a minimum in the energy, giving us zero force
  at that point. Therefore, this interaction needs also a valid ` lj ` file, specified in the input .ffea 
 file through the keyword ` vdw_forcefield_params ` as was explained [earlier](\ref ljPotential).
  

  



