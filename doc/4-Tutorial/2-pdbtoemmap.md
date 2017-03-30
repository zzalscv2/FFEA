
Obtaining an EDM {#pdbtoemmaptut}
=============================

One of the benefits of FFEA is that we can build a system starting from an Electron Density Map, i. e., low resolution EM density data. However, if a higher resolution structure is available then we may wish to begin there instead. From a .pdb file, there are many molecular modelling programs such as VMD, PyMOL and UCSF Chimera with a number of tools for conversion straight to the surface profile that we need. While they are perfectly valid and can be used as input files for the next section of this tutorial, FFEAtools provides it's own set of tools for this purpose guaranteeing that no rotation or translation will affect the resulting map. Thus, in order to generate a pseudo-electron density map, using the FFEAtools type: 

	ffeatools pdbtoemmap [INPUT .pdb fname] [OUTPUT .map fname] [num voxels in x] [num voxels in y] [num voxels in z] [effective atomic radius (angstroms)]

![The atomistic structure of PDB ID:4HEL, visualised using secondary structure.](4hel_2.png "GroEL Atomistic Structure")

The effect of this program is to overlay a cuboidal voxel mesh over the atomistic structure, create a sphere around each atom of the given radius, and assign density to each voxel based on where these spheres (higher density for overlapping spheres and so on). This is an approximation of the surface created by the first electron shell of an atom. [Downloading the PDB file](http://www.rcsb.org/pdb/explore.do?structureId=4HEL) and executing this command in the terminal:

	ffeatools pdbtoemmap 4hel.pdb 4hel_pdbtoem.map 50 50 50 25

Gives the following structure when visualised in Chimera (http://www.cgl.ucsf.edu/chimera/):

![The pseudo-electron density structure of PDB ID:4HEL.](4hel_pdbtoem.png "GroEL Pseudo-Electron Microscopy Structure")

If we compare this to [EMDB ID:EMD-5043](http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-5043), a similar structure, we can see that experimentally obtained electron microscopy data is smoother, easier to work with and therefore more appropriate when just the volumetric data is required and so if low resolution data is available for your protein, use it! 

![The actual electron density structure of EMDB ID:EMD-5043 (recommended contour level 0.42).](emd5043_both.png "GroEL Electron Microscopy Structure")

However, we also see that a small amount of density is contained within the molecule that we will need to get rid of. This is dealt with in the next section.
