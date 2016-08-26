Overview {#overview}
=========

Fluctuating Finite Element Analysis is a new molecular modelling technique, built from the ground-up to support systems that are larger and more complex than those modelled by atomistic molecular dynamics. Instead of modelling biological systems as a collection of connected atoms, it models them as 3D volumes comprised of tetrahedrons. Unlike previous coarse-grained models, the models FFEA generates are viso-elastic continuum solids. Unlike other applications of Finite Element Analysis, these systems are subject to thermal fluctuations.

This technique has the potential to model large, complex systems, made of many molecules, and complex processes at the frontiers of molecular biology. As it does not not require an atomistic level of detail, it can also be used to simulate biological molecules that cannot be imaged using X-ray crystallography.

[cool image]

## Features

 * Protein Interactions:
  * A 6-12 [Lennard-Jones Potential](\ref ljPotential)
  * A repulsive potential that is proportional 
        to the [volume overlap](\ref sPotential).
  * Specific interactions defined using precomputed potentials.
        More documentation can be found [here](\ref fmApproach).
  * Coulombic interactions [EXPERIMENTAL].
 * [Kinetic state changes](\ref kineticApproach) can be simulated together with the continuum model to
    account for conformational changes and binding events.
 * Conversion of EM density data and atomistic structures into FFEA simulations
 * PyMOL visualisation plugin and detailed analysis tools (equilibration, euler characteristic, principal component analysis, gemoetric measurements)


## Videos

[video will go here]

## Publications

   * [A stochastic finite element model for the dynamics of globular macromolecules](http://www.sciencedirect.com/science/article/pii/S0021999112007589),
    R. C. Oliver, D. J. Read, O. G. Harlen, S. A. Harris, J Comp Phys, (2013), 2399:147-165. 
   * List
   * Of
   * Publications

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

## Contribute

FFEA is maintained by a small but dedicated team at the University of Leeds. If you want to see where we can take FFEA, then you can:

   * Use the software for something cool
   * Send bug reports and feature requests to our [issue tracker](https://bitbucket.org/sohpc-ffea/ffea/issues)
   * [Fork us](https://bitbucket.org/sohpc-ffea/ffea/fork)

## FFEA Team

 * Code:
   * Albert Solernou
   * Ben Hanson
   * Robin Richardson
   * [Rob Welch](http://robwel.ch/)
 * Theory
   * [Oliver Harlen](https://www.maths.leeds.ac.uk/index.php?id=263&uid=1025)
   * [Sarah Harris](http://www.comp-bio.physics.leeds.ac.uk/)
   * Robin Oliver
   * [Daniel Read](http://www1.maths.leeds.ac.uk/~djread/)
   * Robin Richardson
   * Ben Hanson
   * Albert Solernou
 * Thanks
   * Lorem Ipsum
   * Dolor Sit
   * Amet



How to read this manual
-----------------------

Because the software is split between the FFEA runner (written in C++) 
  and the rest of FFEA tools (written mostly Python, with some C++ and C), 
  this manual is split also in these two sections: 
   [FFEA_runner](\ref userManual) and [ffeatools](\ref FFEAtools)
  to be used as reference pages.
However, the new user may want start reading the [tutorial](\ref Tutorial),
  and consult the reference pages later. 
