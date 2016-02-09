Basic usage  {#ffea_user}
==========================


Introduction {#ffea_intro}
==========================

In order to run an FFEA trajectory you need a `.vol` mesh for every different 
 `blob` that you want to simulate in your system, and then set up an input 
 file .ffea, where the parameters of the simulations will be described. 

The `.vol` mesh can be generated using 
     [NETGEN](http://sourceforge.net/projects/netgen-mesher/). 
 This mesh will be used as input for an [FFEA_tool](\ref makeffeablob)
 to generate a set of topology files together with an initial .ffea 
 file to be configured by the user. 

Once configured, a trajectory will be run through the command line typing:

    ffea <myInputFile.ffea> 

where ` <myInputFile.ffea> ` has consistently defined all the attributes. 

Finally the trajectory can be analysed using a set of [FFEA_tools](\ref analysisTools).



Input file syntax {#ffea_ifsyntax}
-----------------
In this section we describe the structure of the `.ffea` file. The complete 
 reference to all the keywords describing the system can be found 
 [here](\ref keywordReference).

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

