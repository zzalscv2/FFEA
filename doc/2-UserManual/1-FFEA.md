Basic usage  {#ffea_user}
==========================


Running a trajectory {#ffea_run}
====================


It all starts with:

    ffea <myInputFile.ffea> 

where ` <myInputFile.ffea> ` has consistently defined all the attributes. 



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
  * each ` kinetic ` block defines the structural states and switching rates of each conformation and linear maps between them
  * ` interactions ` is an optional block where interactions are described either in blocks:
      - ` springs ` where springs between nodes are defined or 
      - ` precomp ` where precomputed potentials given in look-up tables.

Attributes within blocks are defined between ` < ` `  > `, using International System Units. 
 For example, kT is written at 300K as:
     
     <kT = 4.11e-21>

Finally, comments are allowed when enclosed between ` &lt;!-- ` and ` --> ` signs, 
  e. g. <!-- this is a comment that does not show off in the HTML version :) --> 
  ` &lt;!--  this is a comment :) --> `.

