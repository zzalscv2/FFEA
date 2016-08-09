Overview {#overview}
=========

[screenshot goes here]

[need a paragraph or so of explanatory text text]

Fluctuating Finite Element Analysis models proteins as visco-elastic bodies subject to
 thermal fluctuations. The model is described in the following paper:

 * [A stochastic finite element model for the dynamics of globular macromolecules](http://www.sciencedirect.com/science/article/pii/S0021999112007589),
    R. C. Oliver, D. J. Read, O. G. Harlen, S. A. Harris, J Comp Phys, (2013), 2399:147-165. 


## Features

 * [Kinetic state changes](\ref kineticApproach) can be simulated together with the continuum model to
    account for conformational changes and binding events.
 * Protein Interactions:
  * A 6-12 [Lennard-Jones Potential](\ref ljPotential)
  * A repulsive potential that is proportional 
        to the [volume overlap](\ref sPotential).
  * Specific interactions defined using precomputed potentials.
        More documentation can be found [here](\ref fmApproach).
  * Coulombic interactions [EXPERIMENTAL].
 * Conversion of EM density data and atomistic structures into FFEA simulations
 * PyMOL visualisation plugin and detailed analysis tools (equilibration, euler characteristic, principal component analysis, gemoetric measurements)


## Videos

## Technology
   * [Boost](http://www.boost.org)
   * [Eigen](http://eigen.tuxfamily.org) (>=3.2.1).   
     FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models.
   * [Doxygen](http://www.doxygen.org) (>= 1.8) [OPTIONAL]   
   * [PyMOL](https://www.pymol.org) (>=1.8) can 
        be used to visualise FFEA systems and trajectories
        as well as molecular and EM systems. Alternatives 
        to visualise molecular systems and create FFEA continuum models
        include [Chimera](https://www.cgl.ucsf.edu/chimera/)
        and [VMD](http://www.ks.uiuc.edu/Research/vmd/).
   * [GTS](http://gts.sourceforge.net) (>=0.7.6)[OPTIONAL]. The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.
   * [NETGEN](https://sourceforge.net/projects/netgen-mesher/) 
   or [TETGEN](http://wias-berlin.de/software/tetgen/) [OPTIONAL]. 
     Programs which convert surface profile into volumetric meshes 
        to be used by FFEA.
   * [pyPcazip](https://pypi.python.org/pypi/pyPcazip) [OPTIONAL]
     Some of the Python FFEA analysis tools interact with these 
     Principal Component Analysis library in order to generate the standard
     PCA output (eigensystems, projections, animations etc)
     obtained from standard from equivalent MD simulations.

How to read this manual
-----------------------

Because the software is split between the FFEA runner (written in C++) 
  and the rest of FFEA tools (written mostly Python, with some C++ and C), 
  this manual is split also in these two sections: 
   [FFEA_runner](\ref userManual) and [FFEA_tools](\ref FFEAtools)
  to be used as reference pages.
However, the new user may want start reading the [tutorial](\ref Tutorial),
  and consult the reference pages later. 
