Surface Profile to Volumetric Mesh {#surftovoltut}
=============================

We now have a completed surface profile, a triangular mesh. Our next aim is to fill this empty surface with tetrahedra, as these base elements form the core of the FFEA algorithm.
Algorithm are available to do this, and our algorithm of choice is provided by NETGEN (http://sourceforge.net/projects/netgen-mesher/), specificallyversion 4.9.13 for stability. However, when passed to the NETGEN volumetric meshing program our current structure gives 448,543 elements. This is almost as bad as having a flully atomistic structure, and so we need to coarsen the surface mesh before creating our volumetric mesh. Luckily, FFEAtools has the appropriate tool for the job:

	FFEA_tools.py surftocgsurf [INPUT .surf fname] [OUTPUT .surf fname] [Length threshold] [Volume conserve (y/n)] [Find smallest edge? (y/n)] OPTIONAL[Coarsening range (xmin, xmax, ymin,ymax, zmin,zmax)]

So we can pass a surface profile, and the program will return a coarsened surface such that the smallest length in the system is greater than the length threshold. Although volume conservation and smallest edge searching are optional, we have found that it is always best to set them active for stability reasons. Also, you can define a box around a sub-section of the structure and coarsening will only occur within that region. So, for example, if we wish to coarsen to a level of 5 angstroms:

	FFEA_tools.py surftocgsurf emd_5043.surf emd_5043_5ang.surf 5 y y

Where we simply use '5' because our original structure was in angstroms. Coarsening to a variety of different length scales gives the following structures:

![Coasened surface profiles of EMDB ID:EMD-5043. From the top, moving clockwise: original, 1A, 3A, 5A, 8A](emd5043_surfcg_ORIG1358.png "GroEL Coarsened Surface Profiles")

A superficial look hints that the 5 angstrom structure may be the best one for the job as it retains all relevent structural detail and hasn't started to become distorted due to the linearisation of it's curvature. However, in general we want to <b> coarsen to the level of the smallest component we are interested in </b>. In the case of GroEL, this may be the central struts, which are still ok in the 8 angstrom structure. So, for computational efficiency, let's go for the 8 angstrom structure!

Once we have our structure, it is a simple matter of asking NETGEN to fill the surface with tetrahedra. This is done as follows:

  * Open NETGEN
  * Import the surface. Click File->Import Mesh... and select a .surf file
  * Fill the structure with tetrahedra. Select 'Generate Mesh' and wait until it finishes
  * Save the volume file. Click File->Save Mesh... and save the file as .vol (we'll call ours emd_5043_8ang.vol)

![Importing a .surf file using NETGEN](netgenimportprintscreen.png "NETGEN - Importing a Surface")
<BR>
![Constructing a .vol file using NETGEN](netgengenmeshprintscreen.png "NETGEN - Building a Volume")
<BR>
![Saving a .vol file using NETGEN](netgensavemeshprintscreen.png "NETGEN - Saving a Volume")

Now we have all of the relevent structural data, all that is left is to convert it into the formats requires for an FFEA simulation!
