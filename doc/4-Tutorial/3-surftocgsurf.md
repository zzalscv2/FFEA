Surface Profile to Coarsened Surface Profile {#surftocgsurftut}
=============================

We now have a completed surface profile, a set of triangles which will form a consistent mesh. But as we saw on the previous page, our model has been generated at quite a high resolution. This is purely a computational effect, designed to exactly preserve the visualisation of the surface we calculated. However, we know that we cannot treat an object with size comparable to individual atoms as a continuum! So we must coarsen our surface to an appropriate level for a continuum mechanical approximation. The GTS libraries (https://sourceforge.net/projects/gts/) can help us with this, as their compilation and installation includes a number of example programs which perform precisely this function. So, we use the following protocol:

	/path/to/gts/binaries/stl2gts < [INPUT .stl fname] > [OUTPUT .gts fname]
	/path/to/gts/source/examples/coarsen -l -c {LENGTH^2} < [INPUT .gts fname] > [OUTPUT .gts fname]
	/path/to/gts/binaries/gts2stl < [INPUT .gts fname] > [OUTPUT .stl fname]

There are many options available to us here, but the above '-l -c {LENGTH^2}' allows us to coarsen based on the length scale of the surface triangles. So, let us use the following values:

	stl2gts < emd_5043.stl > emd_5043.gts
	coarsen -l -c 64 < emd_5043.gts > emd_5043_10ang.gts
	gts2stl < emd_5043_10ang.gts > emd_5043_10ang.stl

Which leaves us with a coarsened surface strucure in which the smallest edge is 10A.


![Surface profile of EMDB ID:EMD-5043 coarsened to 8A resolution.](emd_5043_1angto8angsurf.png "GroEL Coarsened Surface Profile")

As an aside, applications of the coarsen program from GTS to a number of test structures has shown that the volume of the object is conserved by default, which should conserve the elastic strain energy of the object throughout an FFEA simulation i.e. a simulation of an object coarsened to 5A should give approximately the same results as one coarsened to 10A with the same distribution of material parameters. This statement represents a big question in the field of multi-scale modelling and whether it is true in general is under constant research in both our group and others!
   
A superficial look hints that the 5 angstrom structure may be the best one for the job as it retains all relevent structural detail and hasn't started to become distorted due to the linearisation of it's curvature. However, in general we want to <b> coarsen to the level of the smallest component we are interested in </b>. In the case of GroEL, this may be the central struts, which are still ok in the 8 angstrom structure. So, for computational efficiency, let's go for the 8 angstrom structure!
