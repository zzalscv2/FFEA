
Fluctuating Finite Element Method
=================================

Theoretical introduction {#ffea}
========================

Fluctuating Finite Element Analysis models proteins as visco-elastic bodies subject to 
 thermal fluctuations. Details on this approach can be found in this
 [paper](http://www.sciencedirect.com/science/article/pii/S0021999112007589 
         "A stochastic finite element model for the dynamics of globular proteins")

Additionally:

 * Specific interactions can be defined using precomputed potentials. Read the 
    [documentation](\ref fmApproach) if you are interested in doing so.

 * Kinetic state changes can be simulated together with the continuum model to
    account for conformational changes and binding events. Read the 
    [documentation](\ref kineticApproach) if you are interested in doing so.


User documentation  {#ffea_user}
==================

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

