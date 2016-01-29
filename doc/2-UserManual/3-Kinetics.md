
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

 Firstly, within the ` <param> ` block, you need the following extra data:

     <calc_kinetics = 1>
     <kinetics_update = num_steps>
     <kinetics_out_fname = kinetics_output.out>
     <binding_site_params = binding_site.params>

 As well as these, the ` <num_conformations> ` and ` <num_states> ` blocks will now have non-unity values in them.
 
 Secondly, within each ` <blob> ` block, we will want to define multiple ` <conformation> ` blocks corresponding to the
 number defined in the ` <param> ` block for each blob. Optionally, we can now define binding sites within each
 ` <conformation> `:

     <binding_sites = sites.bsites>

 Finally the ` <conformation> ` blocks, we must include the ` <kinetics> ` block. This is as follows:

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

This script requires a completed .ffea input file, and will generate the maps between the structures pointed to by the input file. 
