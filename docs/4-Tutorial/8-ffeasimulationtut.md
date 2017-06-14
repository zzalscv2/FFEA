Running a Simulation {#ffeasimulationtut}
=============================

Before starting the simulation, let's have a look inside the script first.  
It's plain text, so open in your favourite text editor. You should see something like this:

    <param>
       <restart = 0>
       <dt = 1.00e-14>
       <kT = 4.11e-21>
       <check = 10000>
       <num_steps = 1.000000e+11>
       <rng_seed = time>
       <trajectory_out_fname = emd_5043_8ang.ftj>
       <measurement_out_fname = emd_5043_8ang.fm>
       <vdw_forcefield_params = emd_5043_8ang.lj>
       <epsilon = 1.00e-02>
       <max_iterations_cg = 1000>
       <kappa = 1.00e+09>
       <epsilon_0 = 1.00e+00>
       <dielec_ext = 1.00e+00>
       <calc_stokes = 1>
       <stokes_visc = 1.00e-03>
       <calc_vdw = 1>
       <vdw_type = steric>
       <vdw_steric_factor = 0.010000>
       <calc_springs = 0>
       <calc_noise = 1>
       <calc_es = 0>
       <es_update = 1>
       <es_N_x = -1>
       <es_N_y = -1>
       <es_N_z = -1>
       <move_into_box = 1>
       <sticky_wall_xz = 0>
       <wall_x_1 = PBC>
       <wall_x_2 = PBC>
       <wall_y_1 = PBC>
       <wall_y_2 = PBC>
       <wall_z_1 = PBC>
       <wall_z_2 = PBC>
       <es_h = 3>
       <num_blobs = 1>
       <num_conformations = (1)>
    </param>
    
    <system>
       <blob>
          <conformation>
             <motion_state = DYNAMIC>
             <topology = emd_5043_8ang.top>
             <material = emd_5043_8ang.mat>
             <stokes = emd_5043_8ang.stokes>
             <pin = emd_5043_8ang.pin>
             <nodes = emd_5043_8ang.node>
             <surface = emd_5043_8ang.surf>
             <vdw = emd_5043_8ang.vdw>
          </conformation>
          <solver = CG_nomass>
          <scale = 1.00e-10>
       </blob>
    </system>
    
    
The FFEA input file is a structured file that has all the information that FFEA needs 
 to run a simulation. Full detail on the algorithms and options available 
 can be found in the [FFEA runner](\ref userManual) section,
 including the [input file reference](\ref keywordReference), which will also allow
 setting up more sophisticated simulations.
For demonstrating purposes, we need to change a few global parameters. 
 So, for now, go ahead and set `dt = 1e-15`, `num_steps = 1000` and `check = 10`. 
 This will set us up for a run of 1 nanosecond, output a frame every 10 simulation time steps which are 0.01 nanoseconds in length.

Before starting the FFEA runner, you may want to consider setting your OMP_NUM_THREADS environment variable. The recommended number of threads is between 4 and 16. 16 will produce the best results, but will only be around 30-40% faster than running on 4 threads. For now, we will set

	export OMP_NUM_THREADS=4

Now, let's run ffea!

	ffea emd_5043_8ang.ffea

Once the simulation has finished (which should take less than 5 minutes), we can have a look at the structure in the [viewer](\ref FFEAviewertut), and also [perform some analysis](\ref FFEAanalysistut).
A small note; if a simulation crashes with an error talking about elements inverting, not to worry! The simulation time step has been set just a little high. Set ` dt ` to a smaller value, and rerun.


