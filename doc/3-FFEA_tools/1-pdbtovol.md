

FFEA Configuration 1 {#pdbtovol}
=============================

Here we describe how to build and ffea system from experimental and readily available structural data.
For the purposes of this example we will be using the structures PDB ID:4HEL and EMDB ID:EMD-5043, a GroEL protein chaperone.

PDB to EM Density  {#pdbtoemmap}
=============================

One of the benefits of FFEA is that we can build a system starting from low resolution em density data. However, if a higher resolution structure is available then we may wish to begin there instead. From a .pdb file, the Visual Molecular Dynamics (VMD) (http://www.ks.uiuc.edu/Research/vmd/) program contains a number of tools for conversion straight to the surface profile that we need. However, these are often quite tricky to use and so FFEAtools provides it's own set of tools for this purpose. In order to generate a pseudo-electron density map, we use the following program:

	ffeatools pdbtoemmap [INPUT .pdb fname] [INPUT .map fname] [num voxels in x] [num voxels in y] [num voxels in z] [effective atomic radius (angstroms)]

![The atomistic structure of PDB ID:4HEL, visualised using secondary structure.](4hel_2.png "GroEL Atomistic Structure")

The effect of this program is to overlay a cuboidal voxel mesh over the atomistic structure, create a sphere around each atom of the given radius, and assign density to each voxel based on where these spheres (higher density for overlapping spheres and so on). This is an approximation of the surface created by the first electron shell of an atom. Performing this as follows:

	ffeatools pdbtoemmap 4hel.pdb 4hel_pdbtoem.map 50 50 50 25

Gives the following structure when visualised in Chimera (http://www.cgl.ucsf.edu/chimera/):

![The pseudo-electron density structure of PDB ID:4HEL.](4hel_pdbtoem.png "GroEL Pseudo-Electron Microscopy Structure")

If we compare this to EMDB ID:EMD-5043, a similar structure, we can see that experimentally obtained electron microscopy data is smoother, easier to work with and therefore more appropriate when just the volumetric data is required and so if low resolution data is available for your protein, use it! 

![The actual electron density structure of EMDB ID:EMD-5043 (recommended countour level 0.42).](emd5043_both.png "GroEL Electron Microscopy Structure")

However, we also see that a small amount of density is contained within the molecule that we will need to get rid of. This is dealt with in the next section.

EM Density to Surface Profile {#emmaptosurf}
=============================

Once we have a piece of volumetric data such as an em density map, we can create a surface profile of it. However, we first need to make sure that the density map is completely filled, and doesn't contain random offshoots that we don't want (artifacts of the experimental data for example). Our EM density map, as viewed from the top down, containeddensity corresponding to the GroEl molecule being filled when, atomistically, we know that this is not the case. 
