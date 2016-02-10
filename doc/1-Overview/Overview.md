Overview {#overview}
=========


Fluctuating Finite Element Analysis models proteins as visco-elastic bodies subject to
 thermal fluctuations. The model is described in the following paper:

 * [A stochastic finite element model for the dynamics of globular macromolecules](http://www.sciencedirect.com/science/article/pii/S0021999112007589),
    R. C. Oliver, D. J. Read, O. G. Harlen, S. A. Harris, J Comp Phys, (2013), 2399:147-165. 

Additionally:

 * Proteins can interact through:

     - A 6-12 [Lennard-Jones Potential](\ref ljPotential)
     - A repulsive potential that is proportional 
        to the [volume overlap](\ref sPotential).
     - Coulombic interactions [EXPERIMENTAL]. 
     - Specific interactions can be defined using precomputed potentials.
        More documentation can be found [here](\ref fmApproach).

 * Kinetic state changes can be simulated together with the continuum model to
    account for conformational changes and binding events. Read the
    [documentation](\ref kineticApproach) if you are interested in doing so.



General workflow 
----------------


In order to begin an FFEA simulation, you first require a `.vol` mesh for every different
 `blob` that you want to simulate in your system. Once you have the structures, you can set up an input
 file (.ffea) to set the parameters of the simulation and the overall arrangement of your molecules.

The `.vol` mesh can be generated using
     [NETGEN](http://sourceforge.net/projects/netgen-mesher/).
 This mesh is used as input for an [FFEA_tool](\ref makeffeablob)
 to generate a set of structure files together with an initial .ffea
 file to be configured by the user.

Once configured, a trajectory can be simulated through the command line, typing:

    ffea <myInputFile.ffea>

where ` <myInputFile.ffea> ` has consistently defined all the attributes.

Finally the trajectory can be analysed using a set of [FFEA_tools](\ref analysisTools).


How to read this manual
-----------------------

Because the software is split between the FFEA runner (written in C++) 
  and the rest of FFEA tools (written mostly Python, with some C++ and C), 
  this manual is split also in these two sections: 
   [FFEA_runner](\ref userManual) and [FFEA_tools](\ref FFEAtools)
  to be used as reference pages.

However, the new user may want read it following
  the [natural workflow](\ref ffeaWorkFlow):
   * setting up the system with an 
  [FFEA_tool](\ref makeffeablob)
   * configuring the [input file](\ref ffea_iFile) and running an
        FFEA trajectory.
   * and [analysing](\ref analysisTools) the trajectory.



