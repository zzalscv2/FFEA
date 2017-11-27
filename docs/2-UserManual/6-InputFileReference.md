

FFEA Input File Reference {#keywordReference}
=============================================


Param Block {#paramBlock}
=======================
  The param block is mandatory, and 
  starts with the line ` <param> ` and ends with the line ` </param> `. 
In between it can take the following parameters: 


#### System parameters #### 

   * ` restart ` <int> (0) <BR>
        If set up to 1, it will continue from the last snapshot found in 
        ` trajectory_out_fname `. Otherwise, FFEA will start from time-step 0.

   * ` checkpoint_in ` <string> <BR>
        The name of the checkpoint file to be read in case of setting ` <restart = 1> `.
        
   * ` dt ` <float>  <BR>
        Time step. Typically ranging 1e-12 to 1e-14, although the ideal 
        value depends on the shortest edge of the mesh and on the forces 
        acting on it. The higher the forces and the smaller the shortest edge,
        the smaller the time step. There is an 
        [FFEA_tool to calculate the maximum time step](\ref userManual)
        that an mesh can handle. 
   
   * ` num_steps ` <int> <BR>
        Number of time steps to be simulated. 
 
   * ` kT ` <float> <BR>
        The product of the Boltzmann constant k and the temperature T, in ` J `. 

   * ` rng_seed ` <int> <BR>
        The initial seed for the pseudo random number generator. If instead of 
        an integer the string "time" is passed, current time will be used to 
        to get the initial seed.

   * ` epsilon ` <float> (0.01) <BR>
        Error tolerance threshold to determine that Conjugate Gradient has converged.

   * ` max_iterations_cg ` <int> (1000) <BR>
        Maximum number of iterations that the Conjugate Gradient solver will take.
        If the algorithm has not converged, the simulation will be stopped. 
       

   * ` num_blobs ` <int> <BR>
        Number of blobs to be simulated.

   * ` num_conformations ` <(list of ints)> <BR>
        List of comma separated integers enclosed in parenthesis 
        listing the number of conformations
        that exist for each blob. E. g., (2,2,1) would indicate two conformations for
        the first blob, two conformations for the second blob and 1 conformation for 
        the third blob. 

   * ` num_states ` <(list of ints)> <BR>
	List of comma separated integers enclosed in parenthesis defining the number 
        of defined for each blob E. g., (2,2,1) would indicate two conformations for
        the first blob, two conformations for the second blob and 1 conformation for 
        the third blob. A state defines which conformation is active and which binding sites are bound

  

#### Output parameters #### 
  
   * ` check ` <int> <BR>
        FFEA will save coordinates, and energy measurements at every ` check ` steps. 

   * ` trajectory_out_fname ` <string> <BR>
        The name of the file where the coordinates of the trajectory will be saved. 

   * ` measurement_out_fname ` <string> <BR>
        The name of the file where energy measurements will be recorded.

   * ` det_measurement_out_fname ` <string> (Defaulting to ` measurement_out_fname ` with extension replaced with ` .fdm `). <BR>
        The name of the file where detailed measurements will be recorded. 

   * ` checkpoint_out ` <string> (Defaulting to ` ffea-file-name ` with extension replaced with `.fcp`) <BR> 
        The name of the checkpoint file where the details of the state 
          of the last recorded frame will be saved in order to be able to use restarts.

   * ` kinetics_out_fname ` <string> <BR>
        The name of the file where kinetic trajectory will be recorded.

   * ` beads_out_fname ` <string> <BR>
        The name of the file where the trajectory for the beads will be recorded, if found. Currently, restarts are not supported.


#### Enable different calculations #### 

   * ` calc_noise ` <int> (1) <BR>
        Enter 1 or 0 to either enable or disable the effect of the thermal noise.

   * ` calc_stokes ` <int> (1)  <BR>
        1 or 0 

   * ` calc_steric ` <int> (1) <BR>
        Enter 1 or 0 to either enable or disable the [Steric interactions](\ref shortRange).

   * ` calc_ssint ` <int> (0) <BR>
        Enter 1 or 0 to either enable or disable the [Short range interactions](\ref shortRange).

   * ` calc_preComp ` <int> (0) <BR>
        Enter either 1 or 0 to enable or disable the calculation of 
      [pre-computed interactions](\ref fmApproach) between beads.
 
   * ` calc_es `(0) <int> <BR>
        Enter either 1 or 0 to enable or disable the electrostatic interactions. 


#### Steric / Short range forces parameters #### 

   * ` ssint_type ` <string> (steric) <BR>
        Either [lennard-jones](\ref ljPotential), [steric](\ref sPotential), [ljsteric](\ref cPotential) or [gensoft](\ref gPotential) depending on the type of calculations
        to be performed.

   * ` ssint_in_fname ` <string>  <BR>
        The name of the ` .lj ` file that contains the interaction parameters for 
          the short-range interaction potential. Required only if ` ssint_type ` is set to 
          either ` lennard-jones `, ` ljsteric ` or ` gensoft `. More details can be found 
          [here](\ref ljPotential).

   * ` inc_self_ssint ` <int> <BR>
        Enter either 1 or 0 to enable or disable short range interactions 
         between faces within the same blob.

   * ` ssint_cutoff ` see box configuration.

   * ` steric_factor ` <float> (1) <BR>
        Proportionality factor for the steric repulsion approach. More details 
         can be found [here](\ref sPotential).

   * ` steric_dr ` <float> (5e-3) <BR>
        Constant used to calculate the numerial derivative of the steric repulsion.
        The value should only be changed for experienced users trying to use FFEA with 
        elements larger than tens of nm. 

   * ` force_pbc ` <int> (0) <BR>
        Enter either 1 or 0 to enable or disable full periodic boundary conditions specifically for steric, 
        lennard-jones and ljsteric interactions between blobs. Set to 0, wall types set to
        "PBC" will enforce that blobs are translated by a box length in the appropriate direction should 
        their centroid stray out of the box, but does not calculate forces with images of corresponding blobs at 
        the opposite boundary. Set to 1, forces will be correctly calculated with an image of the corresponding 
        blob at the boundaries. Should not be used if blob size greater than half box size in any dimension. Not 
        tested with springs, boundary types other than PBC or electrostatics.

#### Hydrodynamics parameters ####

   * ` stokes_visc ` <float> (1e-3) <BR>
 
#### Electrostatics ####

   * ` es_update `  see box configuration <BR>
        Number of steps after which both the electrostatics energy (if `calc_es`) 
          and the neighbour list are updated.

   * ` dielect_ext ` <float> (1) <BR>
 
   * ` epsilon_0 ` <float> (1) <BR>
        Relative permittivity in the interior of the `blobs`.

   * ` kappa ` <float> (1e9) <BR>
         The inverse of the Debye length in \f$m^{-1}\f$. ` 1e9 ` is a reasonable 
         value for water.

   * ` es_h ` <int> (3) <BR> 
         Number of Debye lengths for which the electrostatic interaction is considered.



#### Simulation box configuration #### 


   * ` es_N_x ` <int> (-1) <BR>
        Number of cells or voxels that the simulation box has in the x direction. If ` es_N_x ` is -1, FFEA will calculate a large enough value to store the system.

   * ` es_N_y ` <int> (-1) <BR>
        Number of cells or voxels that the simulation box has in the y direction. If ` es_N_y ` is -1, FFEA will calculate a large enough value to store the system.

   * ` es_N_z ` <int> (-1) <BR>
        Number of cells or voxels that the simulation box has in the z direction. If ` es_N_z ` is -1, FFEA will calculate a large enough value to store the system.

   * ` ssint_cutoff ` <float> (3e-9) <BR> 
        Length of the side of the voxel.

   * ` es_update ` <int> (10) <BR>
        Number of steps after which the neighbour list is refreshed.


   * ` sticky_wall_xz ` <int> (0) <BR>
         Either 0 or 1, depending on whether the ` xz ` wall, at \f$y=0\f$ is 
         a sticky wall or not. 

   * ` wall_x_1 ` <string> (PBC) <BR>
   * ` wall_x_2 ` <string> (PBC) <BR>
   * ` wall_y_1 ` <string> (PBC) <BR>
   * ` wall_y_2 ` <string> (PBC) <BR>
   * ` wall_z_1 ` <string> (PBC) <BR>
   * ` wall_z_2 ` <string> (PBC) <BR>
     The boundary condition type of each wall. If set to PBC, blobs are translated by a box length in the 
     appropriate direction should their centroid stray out of the box, but does not calculate force interactions 
     with an image of the corresponding blob at the boundary.
        


System Block {#systemBlock}
=======================

  The system block is mandatory starts with the line ` <system> ` and ends with the line ` </system> `. 
It only contains other blocks, without any 'system' parameters.

Blob Block {#blobBlock}
-----------------------

  The blob block starts with the line ` <blob> ` and ends with the line ` </blob> `. At least a `<blob>`
  block is expected. It has both it's own parameters and a set of sub-blocks:


   * ` solver ` <string> (CG_nomass) <BR>
     Specifies the mechanical solver used for all conformations within this block:
	 - **CG**: A solver which builds and inverts the mass matrix of the system, leading to second order Euler integration to update positions (small timestep required)
	 - **CG_nomass**: A solver which builds and inverts the viscosity matrix for a massless system, leading to first order Euler integration to update positions (bigger time-step 		    	     allowed)
	 - **masslumped**: A solver which builds and inverts a purely diagonal mass matrix (hence masslumped) of the system, leading to second order Euler integration to update 			     positions (small time-step required)

   * ` scale ` <float> (1e-10) <BR>
     How much the node positions should be scaled before a simulation begins i.e. if the ` nodes ` file positions are in angstroms, the our ` scale `
     must be 1e-10 to convert everything into meters, the SI units required by FFEA

   * ` centroid ` <list of 3 floats> <BR>
      List of comma separated floats enclosed in parenthesis 
      defining the initial centroid of your blob. 
      This should have the same units as your node file.

   * ` rotation ` either <list of 3 floats> or <list of 9 floats> <BR>
      List of comma separated floats enclosed in parenthesis 
      defining the initial rotation applied to your blob. 
      If you specify 3 values, a rotation is applied first about the x axis, 
      then y, then z.
      If you specify 9 values, this is taken directly as a rotation matrix and applied to your blob.

   * ` velocity ` <list of 3 floats> <BR>
      List of comma separated floats enclosed in parenthesis
      defining the initial velocity of your blob. 
      This has no effect for a system using the CG_nomass ` solver `.

   * ` calc_compress ` <int> (0) <BR>
      Indicates whether to multiply node postions by "compress" variable in all x,y,z during
      initialisation. Does not alter equilibrium dimensions of blob. Not tested with irregularly shaped 
      objects, developed and tested for use with spheres.

   * ` compress ` <float>  <BR>
      Set factor to compress (or expand) blob by. Multiplies all node postions in this blob by this value in all 
      x,y,z during initialisation without changing equilibrium form of blob. Required if calc_compress set to 1.

   * Block `<blob>` expects at least a valid conformation. In the case of multiple conformations for the blob, 
      they have to be defined using multiple `<conformation>` blocks. In the case of a single conformation, 
      the user can choose to use a single `<conformation>` subblock, or to directly include the relevant values 
      in the `<blob>` block. 

### Conformation Block {#conformationBlock} ###

  The conformation block starts with the line ` <conformation> ` and ends with the line ` </conformation> `. 
It contains mostly structural information:

   * ` motion_state ` <enum string> (PBC) <BR>
     How the molecule will move through time:
	DYNAMIC - Fully stochastic motion
	STATIC  - No motion at all (for huge molecules and surfaces)
	FROZEN  - Same as STATIC but all positional data will be printed to trajectories anyway

   * ` nodes ` <string> <BR>
     The file name specifying the node positions

   * ` topology ` <string> [OPTIONAL] <BR>
     The file name specifying the node connectivities (how they form elements). `topology` is only needed if 
     `motion_state` is set to `DYNAMIC`.

   * ` surface ` <string> <BR>
     The file name specifying the surface node connectivities (how they form faces)

   * ` material ` <string> [OPTIONAL] <BR>
     The file name specifying the material properties of each element. 
     ` material ` is only needed if `motion_state` is set to `DYNAMIC`.

   * ` ssint ` <string> [OPTIONAL] <BR>
     The file name specifying the type of ssint interaction for each surface face
     ` ssint ` is only needed if `calc_ssint` is set to 1.

   * ` stokes ` <string> [OPTIONAL] <BR>
     The file name specifying the 'radius' of each node, to calculate local hydrodynamics
     ` material ` is only needed if `motion_state` is set to `DYNAMIC`.

   * ` pin ` <string> <BR>
     The file name specifying which, if any, nodes are pinned in place (to simulate permanent binding, for example)

   * ` binding_sites ` <string> [OPTIONAL] <BR>
     The file name specifying the sets of faces which constitute kinetic binding sites, and the type of site
     ` material ` is only needed if binding sites are being used (read the [kinetics section](\ref kffea_inputfile))

   * ` beads ` <string> [OPTIONAL] <BR>
     The file name specifying the position of every bead in the system. This is only needed on those blobs that 
      use beads for pre-computed potentials, the formatting of this file can be found [here](\ref fm_inputfile).

#### Kinetics Block {#kineticsBlock} ####

  The kinetics block starts with the line ` <kinetics> ` and ends with the line ` </kinetics> `.
It contains information detailing how the different structural conformations kinetically switch between one another: 
 
   * ` states ` <string> <BR>
     The file name specifying the set of states available to the molecule

   * ` rates ` <string> <BR>
     The file name specifying the set of rates at which the molecules will kinetically switch between the defined states

### Maps Block {#mapsBlock} ### 

  The maps block starts with the line ` <maps> ` and ends with the line ` </maps> `.
It contains a list of files detailing how the different structural conformations mathematically relate to one another:

   * ` map (i,j) ` <string> <BR>
     A file name specifying the linear map (non-square matrix) from conformation i to conformation j. The map is highly sparse, and so is 
     stored and applied using the Yale Format for sparse matrices (see [here](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29)).

As a side note, if you have N conformations within your blob, then this section should contain N(N-1) maps (each structure mapping to each other structure, not including itself)

Interactions Block {#interactionsBlock}
---------------------------------------

  The interactions block starts with the line ` <interactions> ` and ends with the line ` </interactions> `.
It only contains other blocks, each of which define a different way blobs can interact between one another.

### PreComp Block {#preCompBlock} ###

  This block is used to describe the pre-computed interactions. It should open with
   the line ` <precomp> ` and close with the line ` </precomp> `. The following 
   fields need to be found in this block if ` calc_preComp ` is set to ` 1 `:

 * ` dist_to_m ` - conversion factor to meters for the distance stored in .pot and .force files.

 * ` E_to_J ` - conversion factor to Joules for the energies stored in .pot files.

 * ` types ` - a comma separated list enclosed in parenthesis with the bead type names, e. g.,
                 (B1, B2, B3).

 * ` inputData ` - can take values 1 and 2 where:
      - 1 will read .force and .pot files
      - 2 will read .pot files and compute the forces.

 The expected format for these files is explained [here](\ref potfile)
 * ` folder ` - path (either absolute or relative to the folder where the .ffea file resides)
                 pointing to the folder storing the .pot (and optionally .force)
                 interaction files.






If the value given for ` types ` were ` (B1, B2, B3) `, FFEA would try to read files:


     B1-B1.pot     B1-B1.force
     B1-B2.pot     B1-B2.force
     B1-B3.pot     B1-B3.force
     B2-B2.pot     B2-B2.force
     B2-B3.pot     B2-B3.force
     B3-B3.pot     B3-B3.force


from the given ` folder `. If a file were not found, the corresponding pair
 would be treated as **inactive**, and a warning would be raised. If no files were
 found, an error would be raised.
 


### Springs Block {#springsBlock} ###

  The springs block starts with the line ` <springs> ` and ends with the line ` </springs> `.
   * ` springs_fname ` <string> <BR>
      The file name containing details of simple linear spring objects between individual nodes in your system.




