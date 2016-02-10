

Input File Reference {#keywordReference}
========================================


Param Block {#paramBlock}
=======================
  The param block starts with the line ` <param> ` and ends with the line ` </param> `. 
In between it can take the following parameters: 


#### System parameters #### 

   * ` restart ` <int> (0) <BR>
        If set up to 1, it will continue from the last snapshot found in 
        ` trajectory_out_fname `. Otherwise, FFEA will start from timestep 0.

   * ` dt ` <float>  <BR>
        Time step. Typically ranging 1e-12 to 1e-14, although the ideal 
        value depends on the shortest edge of the mesh and on the forces 
        acting on it. The higher the forces and the smaller the shortest edge,
        the smaller the time step. There is an 
        [FFEA_tool to calculate the maximum time step](\ref timeStepCalculator)
        that an mesh can handle. 
   
   * ` num_steps ` <int> <BR>
        Number of time steps to be simulated. 
 
   * ` kT ` <float> <BR>
        The product of the Boltzmann constant k and the temprature T, in ` J `. 

   * ` rng_seed ` <int> <BR>
        The initial seed for the pseudo random number generator. If instead of 
        an integer the string "time" is passed, current time will be used to 
        to get the initial seed.

   * ` max_iterations_cg ` <int> (1000) <BR>
        Maximum number of iterations that the Conjugate Gradient solver will take.
        If the algorithm has not converged, the simulation will be stopped. 
       

   * ` num_blobs ` <int> <BR>
        Number of blobs to be simulated.

   * ` num_conformations ` <(list of ints)> <BR>
        List of integers enclosed in parenthesis listing the number of conformations
        that exist for each blob. E. g., (2,2,1) would indicate two conformations for
        the first blob, two conformations for the second blob and 1 conformation for 
        the third blob. 

   * ` num_states ` <(list of ints)> <BR>
	List of integers defining the number of defined for each blob E. g., (2,2,1) would indicate two conformations for
        the first blob, two conformations for the second blob and 1 conformation for 
        the third blob. A state defines which conformation is active and which binding sites are bound

  

#### Output parameters #### 
  
   * ` check ` <int> <BR>
        FFEA will save coordinates, and energy measurements at every ` check ` steps. 

   * ` trajectory_out_fname ` <string> <BR>
        The name of the file where the coordinates of the trajectory will be saved. 

   * ` measurement_out_fname ` <string> <BR>
        The name of the file where energy measurements will be recorded.

   * ` kinetics_out_fname ` <string> <BR>
        The name of the file where kinetic trajectory will be recorded.

#### Enable different calculations #### 

   * ` calc_noise ` <int>  <BR>
        Enter 1 or 0 to either enable or disable the effect of the thermal noise.

   * ` calc_vdw ` <int>  <BR>
        Enter 1 or 0 to either enable or disable the [Short range forces](ref shortRange).

   * ` calc_es ` <int> <BR>
        Enter either 1 or 0 to enable or disable the electrostatic interactions. 

   * ` calc_preComp ` <int> <BR>
        Enter either 1 or 0 to enable or disable the [pre-
 

#### Short range parameters #### 

   * ` vdw_type ` <string> (lennard-jones) <BR>
        Either "lennard-jones" or "steric" depending on the type of calculations
        to be performed.

   * ` vdw_forcefield_params ` <string>  <BR>
        The name of the ` .lj ` file that contains the interaction parameters for 
          the Lennard-Jones potential. More details can be found 
          [here](\ref ljPotential).

   * ` vdw_steric_factor ` <float> <BR>
        Proportionality factor for the steric repulsion approach. More details 
         can be found [here](\ref sPotential).

#### Hydrodynamics parameters ####

   * ` calc_stokes ` <int>  <BR>
        1 or 0 

   * ` stokes_visc ` <float> (1e-3) <BR>
 
#### Electrostatics ####

   * ` es_update ` <int> <BR>
        Number of steps after which the electrostatics energy and the neighbour list
          is refreshed.

   * ` dielect_ext ` <float> (1) <BR>
 
   * ` epsilon ` <float> (0.01) <BR>
        Relative permittivity of the space between ` blobs`.

   * ` epsilon_0 ` <float> (1) <BR>
        Relative permittivity in the interior of the `blobs`.

   * ` kappa ` <float> (1e9) <BR>
         The inverse of the Debye length in \f$m^{-1}\f$. ` 1e9 ` is a reasonable 
         value for water.


#### Simulation box configuration #### 


   * ` es_N_x ` <int> <BR>
        Number of cells or voxels that the simulation box has in the x direction. 

   * ` es_N_y ` <int> <BR>
        Number of cells or voxels that the simulation box has in the y direction. 

   * ` es_N_z ` <int> <BR>
        Number of cells or voxels that the simulation box has in the z direction. 

   * ` es_h ` <int> (3) <BR> 
         The size of a voxel, calculated as ` es_h / kappa` i.e. number of debye lengths.

   * ` sticky_wall_xz ` <int> (0) <BR>
         Either 0 or 1, depending on whether the ` xz ` wall, at \f$y=0\f$ is 
         a sticky wall or not. 

   * ` wall_x_1 ` <enum string> (PBC) <BR>
   * ` wall_x_2 ` <enum string> (PBC) <BR>
   * ` wall_y_1 ` <enum string> (PBC) <BR>
   * ` wall_y_2 ` <enum string> (PBC) <BR>
   * ` wall_z_1 ` <enum string> (PBC) <BR>
   * ` wall_z_2 ` <enum string> (PBC) <BR>

     The boundary condition type of each wall
        


System Block {#systemBlock}
=======================

  The system block starts with the line ` <system> ` and ends with the line ` </system> `. 
It only contains other blocks, and no 'system' parameters.

Blob Block {#blobBlock}
-----------------------

  The blob block starts with the line ` <blob> ` and ends with the line ` </blob> `. 
It has both it's own parameters and a set of subblocks:


   * ` solver ` <enum string> (CG_nomass) <BR>
     Specifies the mechanical solver used for all conformations within this block:
	CG 	   - A solver which builds and inverts the mass matrix of the system, leading to second order Euler integration to update positions (small timestep required)
	CG_nomass  - A solver which builds and inverts the viscosity matrix for a massless system, leading to first order Euler integration to update positions (bigger timestep 		    	     allowed)
	masslumped - A solver which builds and inverts a purely diagonal mass matrix (hence masslumped) of the system, leading to second order Euler integration to update 			     positions (small timestep required)

   * ` scale ` <float> (1e-10) <BR>
     How much the node positions should be scaled before a simulation begins i.e. if the ` nodes ` file positions are in angstroms, the our ` scale `
     must be 1e-10 to convert everthing into meters, the SI units required by FFEA

   * ` centroid ` <list of 3 floats> <BR>
      The initial centroid of your blob. This should have the same units as your node file

   * ` rotation ` <list of 3 floats> / <list of 9 floats> <BR>
      The initial rotation applied to your blob. If you specify 3 values, a rotation is applied first about the x axis, then y, then z.
      If you specify 9 values, this is taken directly as a rotation matrix and applied to your blob.

   * ` velocity ` <list of 3 floats> <BR>
      The initial velocity of your blob. This has no effect for a system using the CG_nomass ` solver `.

### Conformation Block {#conformationBlock} ###

  The conformation block starts with the line ` <conformation> ` and ends with the line ` </conformation> `. 
It contains mostly structural information:

   * ` motion_state ` <enum string> (PBC) <BR>
     How the molecule will move through time:
	DYNAMIC - Fully stochastic motion
	STATIC  - No motion at all (for huge molecules and surfaces)
	FROZEN  - Same as STATIC but all positional data will be printed to trajectories anyway

   * ` nodes ` <string> <BR>
     The filename specifying the node positions

   * ` topology ` <string> <BR>
     The filename specifying the node connectivities (how they form elements)

   * ` surface ` <string> <BR>
     The filename specifying the surface node connectivities (how they form faces)

   * ` material ` <string> <BR>
     The filename specifying the material properties of each element

   * ` vdw ` <string> <BR>
     The filename specifying the type of vdw interaction for each surface face

   * ` stokes ` <string> <BR>
     The filename specifying the 'radius' of each node, to calculate local hydrodynamics

   * ` pin ` <string> <BR>
     The filename specifying which, if any, nodes are pinned in place (to simulate permanent binding, for example)

   * ` binding_sites ` <string> <BR>
     The filename specifying the sets of faces which constitute kinetic binding sites, and the type of site

   * ` beads ` <string> <BR>
     The filename specifying...

#### Kinetics Block {#kineticsBlock} ####

  The kinetics block starts with the line ` <kinetics> ` and ends with the line ` </kinetics> `.
It contains information detailing how the different structural conformations kinetically switch between one another: 
 
   * ` states ` <string> <BR>
     The filename specifying the set of states available to the molecule

   * ` rates ` <string> <BR>
     The filename specifying the set of rates at which the molecules will kinetically switch between the defined states

### Maps Block {#mapsBlock} ### 

  The maps block starts with the line ` <maps> ` and ends with the line ` </maps> `.
It contains a list of files detailing how the different structural conformations mathematically relate to one another:

   * ` map (i,j) ` <string> <BR>

     A filename specifying the linear map (non-square matrix) from conformation i to conformation j. The map is highly sparse, and so is 
     stored and applied using the Yale Format for sparse matrices (https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29)

As a side note, if you have N conforamtions within your blob, then this section should contain N(N-1) maps (each structure mapping to each other structure, not including itself)

Interactions Block {#interactionsBlock}
---------------------------------------

  The interactions block starts with the line ` <interactions> ` and ends with the line ` </interactions> `.
It only contains other blocks, each of which define a different way blobs can interact between one another.

### PreComp Block {#preCompBlock} ###


### Springs Block {#springsBlock} ###

  The springs block starts with the line ` <springs> ` and ends with the line ` </springs> `.
   * ` springs_fname ` <string> <BR>
      The filename containing details of simple linear spring objects between individual nodes in your system.




