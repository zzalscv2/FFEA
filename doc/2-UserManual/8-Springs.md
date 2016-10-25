
Restraints {#springs}
======================

FFEA allows to add restraints between nodes belonging to the same or to different 
 Blobs. These restraints are modelled as springs between the nodes, and so, the user
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
