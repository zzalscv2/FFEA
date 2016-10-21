Volumetric Mesh to FFEA {#voltoffeatut}
=============================

We have a volumetric mesh! This is almost all of the required data needed for an FFEA simulation; all that remains are the material parameters and FFEAtools can help us to set these. The following program enables us to generate a set of FFEA input files that can be immediately used for simulation:

	ffeatools voltoffea --mesh [INPUT .vol fname] --density [INPUT density] --shear-visc [INPUT shear viscosity] --bulk-visc [INPUT bulk viscosity] --shear-mod [INPUT shear modulus] --bulk-mod [INPUT bulk modulus] --stokes-radius [INPUT effective node hydrodynamic radius] --dielectric [INPUT dielectric constant] --make-script

As you can see, this program needs the structural detail of the mesh together with the material properties of our protein in order to generate the file set for FFEA. Let us parametrise GroEL as though it was water, simply for testing purposes (very approximate values obtained from...wikipedia. For testing purposes only!).

Firstly, create a directory called 'simulation' and move into it, and move the .vol file as well:

	mkdir simulation
	cd simulation
	cp ../emd_5043_8ang.vol .

Now, using the .vol file:

	ffeatools voltoffea --mesh emd_5043_8ang.vol --density 1.5e3 --shear-visc 1e-3  --bulk-visc 1e-3 --shear-mod 5.5e8 --bulk-mod 2.2e9 --dielec 1.0 --make-script
	
By default, the program will automatically calculate an approximate value for the total hydrodynamic radius of your object. This can be specified using the `--stokes-radius` argument.

Once the script has run its course, you will see a whole group of new files have appeared:

  * <b>emd_5043_8ang.node</b>		The file containing the positions of each vertex of the mesh
  * <b>emd_5043_8ang.top</b>		The file describing the connectivity of the nodes; how they are arranged into elements
  * <b>emd_5043_8ang.mat</b>		The file holding the material properties of each element in the mesh
  * <b>emd_5043_8ang.surf</b>		The file describing the connectivity of the surface nodes only; how they are arranged into surface faces
  * <b>emd_5043_8ang.stokes</b>		The file containing the effective hydrodynamics radius of each vertex of the mesh
  * <b>emd_5043_8ang.vdw</b>		The file containing the van der Waals type of each face (more on this later)
  * <b>emd_5043_8ang.pin</b>		The file containing a simple list of all of the 'pinned' nodes, i.e. nodes which will not move during the simulation
  * <b>emd_5043_8ang.lj</b>		The file containing the various parameter sets for van der Waals interactions between active faces (again, more on this later) 
  * <b>emd_5043_8ang.ffea</b>		The script file used by FFEA which points to all of the above files!

Full description of these files can be found [here](\ref ioFiles). Now, emd_5043_8ang.ffea is ready to run immediately!

