Coarsening a Surface Profile {#surftocgsurftut}
=============================

We now have a completed surface profile, a set of triangles which will form a consistent mesh. But as we saw on the previous page, our model has been generated at quite a high resolution. This is purely a computational effect, designed to exactly preserve the visualisation of the surface we calculated. However, we know that we cannot treat an object with size comparable to individual atoms as a continuum! So we must coarsen our surface to an appropriate level for a continuum mechanical approximation. Although the GTS libraries are the standard for this, FFEA has its own set of coarsening tools that produce more reliable models.

The FFEA surface coarse grainer is used like this:

	ffeatools surftocgsurf [Input filename] [Output filename] [Coarseness level (A)] [Volume conserve (y/n)] [Find smallest edge? (y/n)] OPTIONAL[Coarsening range (xmin, xmax, ymin,ymax, zmin,zmax)]

The recommended settings are to have volume conserve and find smallest edge turned on. The coarseness level is discussed below.

To coarsen the surface, run

	ffeatools surftocgsurf emd_5043.obj emd_5043_8ang.stl 8 y y

Which leaves us with a coarsened surface strucure in which the smallest edge is 8A. We also convert it to the .stl format.


![Surface profile of EMDB ID:EMD-5043 coarsened to 8A resolution.](emd_5043_1angto8angsurf.png "GroEL Coarsened Surface Profile")

> As an aside, applications of the program ` coarsen ` (part of the examples of the GTS library) to a number of test structures has shown that the volume of the object is conserved by default, so it should be a valid alternative to our ffeatools.  Volume conservation is important as it is related to the elastic strain energy of the object throughout an FFEA simulation i.e. a simulation of an object coarsened to 5A should give approximately the same results as one coarsened to 8A with the same distribution of material parameters. This statement represents a big question in the field of multi-scale modelling and whether it is true in general is under constant research in both our group and others!
   

A superficial look hints that a 5 angstrom structure may be the best one for the job as it retains all relevent structural detail and hasn't started to become distorted due to the linearisation of its curvature. However, in general we want to **coarsen to the level of the smallest component we are interested in**. In the case of GroEL, this may be the central struts, which are still ok in the 8 angstrom structure. So, for computational efficiency, let's go for the 8 angstrom structure!
