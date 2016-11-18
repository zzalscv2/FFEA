Restraints {#restraints}
======================

FFEA allows to add different types of restraints, through [springs](\ref springs), 
 and [pinned nodes](\ref pinnednodes) and [frozen blobs](\ref frozen). The following 
 sections go through the detail on what is computed in every case.
 

Springs {#springs}
==================

Adding distance restraints between nodes belonging to the same or to different 
 Blobs is done in FFEA using springs between the nodes. In that case, the user
 is able to specify the strength and equilibrium distance of the Hookean potential. 
In order to use restraints, one has set:

    < calc_springs = 1 > 

in the input .ffea file. In addition to that, ` <system> ` has to include the 
 ` <interactions> ` block with the ` <springs> ` block, where the ` springs_fname ` 
 is specified. Thus, a valid input ffea file simulating two blobs with 
 restraints would have a structure:

     <param>
        ...
        <calc_springs=1>
        ...
     </param>
     <system>
        <blob>
           ...
        </blob>
        <blob>
           ...
        </blob>
        <interactions>
           <springs>
              <springs_fname=myspringsfile.springs>
           </springs>
        <interactions> 
     </system> 
 

It is mandatory that the ` springs_fname ` suffix is "springs". Full description
 of the format of the input ` springs_fname ` can be found [here](\ref ifsprings).


Freezing degrees of freedom {#freezeing}
===========================
Degrees of freedom can be individually or collectively removed through pinned nodes 
 or freezing blobs. In any case, the existing active faces within the frozen structures
 as well as those beads assigned to the corresponding elements will remain involved 
 in the calculation of forces that will be applied to the rest of the structures.

Pinned nodes  {#pinnednodes}
-----------
If using a list of pinned nodes through [pin](\ref conformationBlock) one has to pass
 in a [.pin file](\ref ifpin) with a list of nodes. The position of the nodes within
 that list will never change during simulation.  


Frozen blob {#frozen}
-----------
One can freeze a whole blob passing the flag ` STATIC ` or ` FROZEN ` to the ` motion_state ` 
 of a Blob [conformation](\ref conformationBlock). Again, the position of any of the nodes 
 within this conformation will never be updated.




