

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


Blob Block {#blobBlock}
-----------------------


### Conformation Block {#conformationBlock} ###


#### Kinetics Block {#kineticsBlock} ####


### Maps Block {#mapsBlock} ### 


Interactions Block {#interactionsBlock}
---------------------------------------


### PreComp Block {#preCompBlock} ###


### Springs Block {#springsBlock} ###





