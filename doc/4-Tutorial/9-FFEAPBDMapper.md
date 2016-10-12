

FFEA to PDB Mapping {#FFEAPDBMapper}
===============================

With `multiscale' being somewhat of a buzzword at the moment in the molecular modelling community, albeit a valid one,
we have always had the underlying atomic structure of our continuum models in mind as we develop our simulations. As scuh,
we have developed a procedure to transform from the continuum structure back to the underlying atomic structure so that following
an FFEA simulation, further atomistic simulations can be performed.

We emphaisise that this is only possible if the atomic structure of your model has been experimentally calculated in the first place.
We are not magic...

Mapping Theory {#FFEApdbmappertheory}
===============================

If we build our FFEA structure from the underlying PDB structure, then at the end of the initialisation routines the continuum and atomistic
structures should be aligned with one another (unless some global translation / rotation has been applied, which can be reversed). When spatially overlapped
like this, the constant topology of an FFEA structure allows us to generate a linear mapping between the node of the mesh and the atoms of the PDB, which can
be written as a matrix,

	\f$\vec{x}^{atoms} = \mathbf{M}\vec{x}^{nodes}\f$

Where \f$\mathbf{M}\f$ is a non-square matrix. We can constuct this matrix by making the position of every atom, \f$\vec{x}_{\alpha}^{atoms}\f$, a function of only the nodes forming it's containing element, or nearest element,

	\f$\vec{x}^{atoms}_{\alpha} = m_{\alpha \gamma} \vec{x}^{nodes}_{\gamma}\f$

where \f$m_{\alpha \gamma}\f$ is a smaller, \f$4 \times 4\f$ sub-matrix of \f$\mathbf{M}\f$. Using the shape functions of the element formed by the nodes \f$\gamma \f$, we can find the components of \f$m_{\alpha \gamma}\f$ and use it to populate \f$\mathbf{M}\f$. This matrix, which is based on the equilibrium structures of both the FFEA and PDB models, can be used to convert from a non-equilibrium FFEA structure (following a simulation) to the approximately equivalent atomistic structure.

Mapping Practicality {#FFEApdbmapperpractical}
===============================  

The FFEAtools package contains a number of scripts to allow the creation and application of this map. For an initialised FFEA system and its associated atomistic PDB, we use the front end of the toolkit to access the script:

	
	python /path/to/ffea/lib/pythonX.Y/ffeatools.py makestructuremap -i FFEAstructure.node -t FFEAstructure.top -o Atomisticstructure.pdb -m FFEAtoatoms.map

This will use the [above method](\ref #FFEApdbmappertheory) to generate your mapping script. This is a highly sparse matrix, so we also require a conversion to the Yale sparse matrix format to save space. This is currently implemented as another script (sorry!):

	
	python /path/to/ffea/lib/pythonX.Y/ffeatools.py maptosparse FFEAtoatoms.map FFEAtoatoms_sparse.map

And there we go, we have the map. This matrix can be applied to any simulation frame calculated by FFEA in order to generate an atomistic structure. FFEA tools again provides a script for this procedure:

	python /path/to/ffea/lib/pythonX.Y/ffeatools.py maptraj FFEAtrajectory.ftj FFEAtoatomstrajectory.pdb FFEAtoatoms_sparse.map Atomisticstructure.pdb

The above script converts the entire FFEA trajectory into an atomistic one. This is done by applying the map to the FFEA trajectory (which has the same node order as the original .node file), then importing the original .pdb topology onto the new atomic positions (again, same order as the original) and writing the whole thing to a file.

Visualisation in Pymol of both this new atomistic trajectory and the FFEA trajectory will show clearly how the mapping procedure work, but it must be emphasised that this is an entirely non-physical tranformation. If atomic positions are extrapolated by the process, then bond lengths may be much greater than they should be. Interpolation of atomic positions based on the FFEA simulation can cause errors in the dihedral angles and other angular properties. If you wish to perform atomistic simulations on the back of this mapping procedure, we advise that you minimise and equilibrate the resulting atomistic structure as best you can first. 


