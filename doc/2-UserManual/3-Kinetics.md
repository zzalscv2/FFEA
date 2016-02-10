
Kinetic FFEA {#kineticApproach}
======================

Theoretical introduction {#kffea}
========================

One of the main advantages of FFEA compared with more common techniques such as Molecular Dynamics (MD)
 is the larger simulation time it is able to reach. At the &micro;s regime and beyond, certain proteins,
 such as molecular motors, can undergo internal conformational changes and external binding events. These
 processes are themselves atomic in nature, triggered by nanoscale processes such as ATP hydrolysis, yet 
 the net result occurs at the mesoscale and so are out of range of MD. We have enabled the simulation of these
 events by integrating FFEA with a kinetic model, enabling switching between different conformational
 states and bound / unbound states whilst simulataneously continuing the stochastic dynamic simulation.

Relevant fields at the input file {#kffea_inputfile}
=================================
In order to implement kinetics, the following extra data is required by the simulation.

 Firstly, we need to tell FFEA that we want to run with kinetics. Within the ` <param> ` block, we need the following extra parameters:

 To initialise kinetics:
     <calc_kinetics = 1>

 How often we want to update our kinetic states:
     <kinetics_update = num_steps>

 The name of the output file containing kinetic data:
     <kinetics_out_fname = kinetics_output.out>

 And that file defining how binding sites interact (if they exist):
     <binding_site_params = binding_site.params>

 As well as these, the ` <num_conformations> ` and ` <num_states> ` blocks will now have non-unity values in them.
 
 We now need to define our kinetic states. Within each ` <blob> ` block, we may want to define multiple ` <conformation> ` blocks corresponding to the
 number defined in the ` <param> ` block for each blob. These conformations each contain independent structural data which we kinetically transform between.

 Additionally, we can define binding sites within each ` <conformation> `:

     <binding_sites = sites.bsites>

 Finally the ` <blob> ` blocks must include the ` <kinetics> ` block as well as the ` <conformation> ` block, to define how often each conformation maps onto another. 
 This is as follows:

     <kinetics>
         <states = states_fname.states>
         <rates = rates_fname.rates>
         <maps>
           ...
          </maps>
     </kinetics>

 The ` <states = ...> ` tag defines what conformations are active and what binding sites are bound together
 in each state.  The ` <rates = ...> ` tag defines the rates at which these states switch between one another.

 For every pair of conformations, two structural maps are required in order to switch between the pairs of structures:

     <maps>
         <map (conf_index_A, conf_index_B) = AtoB.map>
         <map (conf_index_B, conf_index_A) = BtoA.map>
     </maps>

Implementation details {#kffea_implementation}
----------------------
In order to generate the maps between structures, a script is provided within the FFEA toolkit:

     FFEAtools makekineticmaps

This script requires a completed and stable .ffea input file, and will generate the maps between any two structures pointed to by the input file. 
The algorithm requires continuous user input the whole way through due to the not-standard method of defining maps (i.e. maximising volume overlap may not be valid).

