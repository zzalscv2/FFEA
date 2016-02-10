
FFEA Input File  {#ffea_iFile}
==============================

In this section we describe the structure of the `.ffea` file. The complete 
 reference to all the keywords describing the system can be found 
 [here](\ref keywordReference). More information on every simulation approach
 will be found in the corresponding pages.



Input file syntax {#ffea_ifsyntax}
-----------------

The input file has two main blocks: 
  * param: describing global parameters of the system
  * system: describing the system itself.

Blocks open with: ` <block> ` and closed with ` </block> `. Therefore a valid input file
 will look like:


    <param>
      ...
    </param>

    <system>
        <blob>
            <conformation>
             ...
            </conformation>

	    <kinetics>
	     ...
	    </kinetics>
        </blob>

        <interactions> 
           ...
         </interactions>

    </system>     

where:
  * each ` blob ` block defines a protein in the simulation
  * each ` conformation ` describes structurally a stable conformation of a protein
  * each ` kinetic ` block defines the how multiple structural states can kinetically switch between one another
  * ` interactions ` is an optional block where interactions between blobs are defined, containing at least one of the following blocks:
      - ` springs ` where springs between nodes are defined or 
      - ` precomp ` where precomputed potentials are given in look-up tables.

Attributes within blocks are defined between ` < ` `  > `, using International System Units. 
 For example, temperature is set to 300K using kT, as:
     
     <kT = 4.11e-21>

 The simulation length is set using num_steps and will be cast to an integer:

     <num_steps = 1e8>

 Filenames can be written relative to the current working directory or absolute:

    <trajectory_out_fname = /path/to/your/file/trajectory.out>

Finally, comments are allowed when enclosed between ` &lt;!-- ` and ` --> ` signs, 
  e. g. <!-- this is a comment that does not show off in the HTML version :) --> 
  ` &lt;!--  this is a comment :) --> `.

