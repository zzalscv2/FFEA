
FFEA Input File  {#ffea_iFile}
==============================

In this section we describe the structure of the `.ffea` file. The complete 
 reference to all the keywords describing the system can be found 
 [here](\ref keywordReference). More information on every simulation approach
 will be found in the corresponding pages.



Input file syntax {#ffea_ifsyntax}
==================================

The input file has two main blocks: 
  * param: describing global parameters of the system
  * system: describing the system itself.

Blocks open with: ` <block> ` and closed with ` </block> `, and the structure of a minimal
 valid ` .ffea ` file looks like:


    <param>
      ...
    </param>

    <system>
        <blob>
         ... 
        </blob>
    </system>     

where:
  * the `<param>` block contains parameters affecting the whole `system`.
  * the `<system>` block contain subblocks describing the blobs and possibly their interactions.
  * each ` blob ` block defines a protein in the simulation

Attributes within blocks are defined between ` < ` `  > `, using International System Units. 
 For example, temperature is set to 300K using kT, as:
     
     <kT = 4.11e-21>

 The simulation length is set using num_steps and will be cast to an integer:

     <num_steps = 1e8>

 File names can be written relative to the directory where the .ffea script resides, 
or giving an absolute path:

    <trajectory_out_fname = /path/to/your/file/trajectory.out>

Comments are allowed when enclosed between ` &lt;!-- ` and ` --> ` signs, 
  e. g. <!-- this is a comment that does not show off in the HTML version :) --> 
  ` &lt;!--  this is a comment :) --> `.

A minimal input file is generated automatically through [voltoffea](\ref voltoffeatut),
 the next sections describe how to set up more complex systems, and all the 
 keywords and blocks are listed and described in the [reference site](\ref keywordReference). 

 
 
