

FFEA Configuration 1 {#pdbtovol}
=============================

Here we describe how to build and ffea system from experimental and readily available structural data.
For the purposes of this example we will be using the structures PDB ID:4HEL and EMDB ID:EMD-5043, a GroEL protein chaperone.

PDB to EM Density  {#pdbtovol}
=============================

One of the benefits of FFEA is that we can build a system starting from low resolution em density data. However, if a higher resolution structure is available then we may wish to begin there instead. From a .pdb file, the Visual Molecular Dynamics (VMD) program contains a number of tools for conversion straight to the surface profile that we need. However, these are often quite tricky to use and so FFEAtools provides it's own set of tools for this purpose. In order to generate a pseudo-electron density map, we use the following program:

	FFEA_tools.py pdbtoemmap [INPUT .pdb fname] [INPUT .map fname] [num voxels in x] [num voxels in y] [num voxels in z] [effective atomic radius (angstroms)]


