Volumetric Mesh to FFEA {#voltoffeatut}
=============================

We have a volumetric mesh! This is almost all of the required data needed for an FFEA simulation; all that remains are the material paramters and FFEAtools can help us to set these. The following program enables us to generate a set of FFEA input files that can be immediately used for simulation:

	ffeatools voltoffea -mesh [INPUT .vol fname] -density [INPUT density] -shear_visc [INPUT shear viscosity] -bulk_visc [INPUT bulk viscosity] -shear_mod [INPUT shear modulus] -bulk_mod [INPUT bulk modulus] -stokes_radius [INPUT effective node hydrodynamic radius] -dielec [INPUT dielectric constant] -make_script

As you can see, this program needs the structural detail of the mesh together with the material properties of our protein in order to generate the file set for FFEA. Let us parameterise GroEL as though it was water, simply for testing purposes (very approximate values obtained from...wikipedia. For testing purpses only!).

Firstly, create a directory called 'simulation' and move into it, and move the .vol file as well:

	mkdir simulation
	cd simulation
	cp ../emd_5043_10ang.vol .

Now, using the .vol file:

	ffeatools voltoffea -mesh emd_5043_10ang.vol -density 1.5e3 -shear_visc 1e-3  -bulk_visc 1e-3 -shear_mod 5.5e8 -bulk_mod 2.2e9 -dielec 1.0 -make_script

The program will ask if you want a stokes_radius calculating. Select 'y' and the program calculates an approximate value for the total hydrodynamic radius of your object. Now you will see a whole group of new files have appeared:

  * <b>emd_5043_10ang.node</b>		The file containing the positions of each vertex of the mesh
  * <b>emd_5043_10ang.top</b>		The file describing the connectivity of the nodes; how they are arranged into elements
  * <b>emd_5043_10ang.mat</b>		The file holding the material properties of each element in the mesh
  * <b>emd_5043_10ang.surf</b>		The file describing the connectivity of the surface nodes only; how they are arranged into surface faces
  * <b>emd_5043_10ang.stokes</b>		The file containing the effective hydrodynamics radius of each vertex of the mesh
  * <b>emd_5043_10ang.vdw</b>		The file containing the Van der Waals type of each face (more on this later)
  * <b>emd_5043_10ang.pin</b>		The file containing a simple list of all of the 'pinned' nodes, i.e. nodes which will not move during the simulation
  * <b>emd_5043_10ang.lj</b>		The file containing the various parameter sets for Van der Waals interactions between active faces (again, more on this later) 
  * <b>emd_5043_10ang.ffea</b>		The script file used by FFEA which points to all of the above files!

Now, emd_5043_10ang.ffea is ready to run immediately! But let's have a look inside the script first. It's plain text, so open in your favourite text editor. You should see something like this:

![An example of an FFEA script file (.ffea)](ffeascript.png "GroEL FFEA Script")

The structure of file script tries to be as helpful as possible, but if you are unsure what is going on, check the [input file reference](\ref keywordReference) for details. What we need to do first of all is change a few global parameters. So, go ahead and set 'dt = 1e-13', 'num_steps = 10000' and 'check = 100'. This will set us up for a run of 1 nanosecond, output a frame every 100 simulation timesteps which are 0.01 nanoseconds in length. Now, let's run ffea!

	ffea emd_5043_10ang.ffea

Once the simulation has finished (which should take about 5 minutes), we can have a look at the structure in the [viewer](\ref FFEAviewertut), and also [perform some analysis](\ref FFEAanalysistut).
A small note; if the simulation crashes with an error talking about elements inverting, not to worry! Our simulation timestep has been set just a little high. Set dt = 5e-13/1e-14 or even smaller if you wish, and rerun.




