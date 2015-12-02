
Fluctuating Finite Element Method
=================================

Theoretical introduction {#ffea}
========================

Fluctuating Finite Element Analysis models proteins as visco-elastic bodies subjected to 
 thermal fluctuations. Details on this approach can be found in this
 [paper](http://www.sciencedirect.com/science/article/pii/S0021999112007589 
         "A stochastic finite element model for the dynamics of globular proteins")

Specific interactions can be defined using precomputed potentials. Read the 
 [documentation](\ref fmApproach) if you are interested in using so.



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
        </blob>

        <interactions> 
           ...
         </interactions>

    </system>     

where:
  * each ` blob ` block define a protein in the simulation
  * each ` conformation ` describes structurally a stable conformation of a protein
  * ` interactions ` is an optional block where interactions are described either in blocks:
      - ` springs ` where springs between nodes are defined or 
      - ` precomp ` where precomputed potentials given in look-up tables.

Attributes within blocks are defined between ` < ` `  > `, using International System Units. 
 For example, KT is written at 300K as:
     
     <KT = 4.11e-21>

Finally, comments are allowed when enclosed between ` &lt;!-- ` and ` --> ` signs, 
  e. g. <!-- this is a comment that does not show off in the HTML version :) --> 
  ` &lt;!--  this is a comment :) --> `.

