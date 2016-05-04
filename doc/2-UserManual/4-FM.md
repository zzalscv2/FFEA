
Precomputed potentials {#fmApproach}
======================


Relevant fields at the input file {#fm_inputfile}
=================================
In order to use pre-computed potentials one needs to enter data in three different blocks.
 Firstly, you need to state:

     <calc_preComp = 1> 

within the ` <param> ` block. 

Secondly, beads for each conformation 
 need to be in a separate file using the standard PDB format:


    ATOM      1   B4 Pro A   1    -142.098 -12.009   5.325
    ATOM      2   B1 Pro A   1    -145.067  -7.988   7.025
    ATOM      3   B1 Pro A   1    -146.489 -11.901   3.324
    ATOM      4   B4 Pro A   1    -151.897  -7.813   7.042
    ATOM      5   B2 Pro A   1    -154.550 -11.857   6.645
    ATOM      6   B2 Pro A   1    -155.957  -6.345  10.220
    ATOM      7   B4 Pro A   1    -157.703  -7.778   6.177


where the only fields that are taken into account are 
   ` B1 `, ` B2 `, ... standing for the bead type (single word, no further restrictions)
 and the position of each bead ` x `, `y`, `z` 
 (using the same units as the ` nodes ` file). 
 This file is entered as ` beads ` within block ` <conformation> `:

     <beads = beadsForBlobXConfY.pdb> 

At the begining of the simulation, each bead is assigned onto a tetrahedron, and 
 during the simulation, the forces that these beads experiment will be linearly
 interpolated onto the nodes of the corresponding tetrahedra. The initial assignment
 is done to the nearest tetrahedra centroid. However, this can be altered individually
 adding a directive at the end of every line within the "beads" file. More specifically,
 if we were assigning bead 2, to a tetrahedron that has either nodes 7, 10 or 15 to 20, 
 we would write:


    ATOM      1   B4 Pro A   1    -142.098 -12.009   5.325
    ATOM      2   B1 Pro A   1    -145.067  -7.988   7.025 <nodes = 7,10,15-20>
    ATOM      3   B1 Pro A   1    -146.489 -11.901   3.324
    ATOM      4   B4 Pro A   1    -151.897  -7.813   7.042
    ATOM      5   B2 Pro A   1    -154.550 -11.857   6.645
    ATOM      6   B2 Pro A   1    -155.957  -6.345  10.220
    ATOM      7   B4 Pro A   1    -157.703  -7.778   6.177


where the directive ` nodes ` ensures that no accidental confusion arises with
 B-factors or other information after the position of the beads. 


Finally, the rest of the input flags are passed within block ` <precomp> ` in ` <interactions> `. 
 Fields to be provided are: 
 * ` types ` - a comma separated list between parenthesis with the bead type names, e. g., 
                 (B1, B2, B3, B4).
 * ` inputData ` - can take values 1 and 2 where:
      - 1 will read .force and .pot files 
      - 2 will read .pot files and compute the forces. 

     The expected format for these files is explained [here](\ref potfile)

 * ` folder ` - relative or absolute path to the folder storing the .pot 
                 (and optionally .force) interaction files. 
 * ` approach ` - must be ` solid `. This field will be removed.
 * ` dist_to_m ` - conversion factor to meters for the distance stored in .pot and .force files. 
 * ` E_to_J ` - conversion factor to Joules for the energies stored in .pot and .force files. 


Implementation details {#fm_implementation}
----------------------

 I try to keep all the methods and functions together in a module named ForceMatch.




