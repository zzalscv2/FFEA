Volumetric Mesh to FFEA {#voltoffeatut}
=============================

We have a volumetric mesh! This is almost all of the required data needed for an FFEA simulation; all that remains are the material paramters and FFEAtools can help us to set these. The following program enables us to generate a set of FFEA input files that can be immediately used for simulation:

	ffeatools voltoffea --mesh [INPUT .vol fname] --density [INPUT density] --shear_viscosity [INPUT shear viscosity] --bulk_viscosity [INPUT bulk viscosity] --shear_modulus [INPUT shear modulus] --bulk_modulus [INPUT bulk modulus] --stokes_radius [INPUT effective node hydrodynamic radius] --dielectric [INPUT dielectric constant] --make_script

As you can see, this program needs the structural detail of the mesh together with the material properties of our protein in order to generate the file set for FFEA. Let us parameterise GroEL as though it was water, simply for testing purposes (very approximate values obtained from...wikipedia. For testing purpses only!).

Firstly, create a directory called 'simulation' and move into it, and move the .vol file as well:

	mkdir simulation
	cd simulation
	cp ../emd_5043_8ang.vol .

Now, using the .vol file:

	ffeatools voltoffea --mesh emd_5043_8ang.vol --density 1.5e3 --shear_visc 1e-3  --bulk_visc 1e-3 --shear_mod 5.5e8 --bulk_mod 2.2e9 --dielec 1.0 --make_script
	
By default, the program will automatically calculate an approximate value for the total hydronamic radius of your object. This can be specified using the `--stokes_radius` argument.

Once the script has run its course, you will see a whole group of new files have appeared:

  * <b>emd_5043_8ang.node</b>		The file containing the positions of each vertex of the mesh
  * <b>emd_5043_8ang.top</b>		The file describing the connectivity of the nodes; how they are arranged into elements
  * <b>emd_5043_8ang.mat</b>		The file holding the material properties of each element in the mesh
  * <b>emd_5043_8ang.surf</b>		The file describing the connectivity of the surface nodes only; how they are arranged into surface faces
  * <b>emd_5043_8ang.stokes</b>		The file containing the effective hydrodynamics radius of each vertex of the mesh
  * <b>emd_5043_8ang.vdw</b>		The file containing the Van der Waals type of each face (more on this later)
  * <b>emd_5043_8ang.pin</b>		The file containing a simple list of all of the 'pinned' nodes, i.e. nodes which will not move during the simulation
  * <b>emd_5043_8ang.lj</b>		The file containing the various parameter sets for Van der Waals interactions between active faces (again, more on this later) 
  * <b>emd_5043_8ang.ffea</b>		The script file used by FFEA which points to all of the above files!

Now, emd_5043_8ang.ffea is ready to run immediately! But let's have a look inside the script first. It's plain text, so open in your favourite text editor. You should see something like this:

![An example of an FFEA script file (.ffea)](ffeascript.png "GroEL FFEA Script")

The structure of file script tries to be as helpful as possible, but if you are unsure what is going on, check the [input file reference](\ref keywordReference) for details. What we need to do first of all is change a few global parameters. So, for now, go ahead and set 'dt = 1e-15', 'num_steps = 1000' and 'check = 10'. This will set us up for a run of 1 nanosecond, output a frame every 10 simulation timesteps which are 0.01 nanoseconds in length.

Before starting the FFEA runner, you may want to consider setting your OMP_NUM_THREADS environment variable. The recommended number of threads is between 4 and 16. 16 will produce the best results, but will only be around 30-40% faster than running on 4 threads. For now, we will set

```sh
export OMP_NUM_THREADS=4
```

Now, let's run ffea!

	ffea emd_5043_8ang.ffea

Once the simulation has finished (which should take about 5 minutes), we can have a look at the structure in the [viewer](\ref FFEAviewertut), and also [perform some analysis](\ref FFEAanalysistut).
A small note; if the simulation crashes with an error talking about elements inverting, not to worry! Our simulation timestep has been set just a little high. Set dt = 5e-13/1e-14 or even smaller if you wish, and rerun.




