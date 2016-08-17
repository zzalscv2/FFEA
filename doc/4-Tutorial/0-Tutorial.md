Tutorial {#Tutorial}
=========================

Overview 
----------------

In order to begin an FFEA simulation, you first require a `.vol` mesh for every different
 `blob` that you want to simulate in your system. Once you have the structures, you can then set up an input
 file (.ffea) to set the parameters of the simulation and the overall arrangement of your molecules.

The `.vol` mesh can be generated using
     [NETGEN](http://sourceforge.net/projects/netgen-mesher/).
 This mesh is used as input for an [FFEA_tool](\ref makeffeablob)
 to generate a set of structure files together with an initial .ffea
 file to be configured by the user.

Once configured, a trajectory can be simulated through the command line, typing:

    ffea <myInputFile.ffea>

where ` <myInputFile.ffea> ` has consistently defined all the attributes.

Finally the trajectory can be analysed using a set of [ffeatools](\ref analysisTools). More details on all these steps can be found in the following chapters
 of this manual.

Content:

- @subpage pdbtoemmaptut
- @subpage emmaptosurftut
- @subpage surftocgsurftut
- @subpage surftovoltut
- @subpage voltoffeatut
- @subpage FFEAviewertut
- @subpage ffeasimulationtut
- @subpage FFEAanalysistut

