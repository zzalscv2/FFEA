FFEA runner {#userManual}
=========================

The FFEA runner is the program that will compute a FFEA trajectory given 
 an initial system under a set of conditions, all contained in the
 ` FFEA Input File `. The program is meant to run from the command line 
 typing:
 
      ffea <myInputFile.ffea>


A number of optional arguments are available:

	'-h/--help' - Prints out a help text to the command prompt
	'-v/--verbose' - Prints out to the command line a user defined level of information
	'-d/--no-detail' - Only writes measurement data to files for the total global system, not each individual blob and pair of blobs (the .fdm file is not produced)
	'-m/--mode' - Select the mode of FFEA you would like to run: Full simulation [default] (0), Linear Elastic Model (1), Linear Dynamic Model (2) or Timestep Calculator (3). 
 
FFEA can run in parallel using multiple threads through OpenMP. By default, 
 it will try to use as many threads as cores found. One can control the 
 number this benaviour adjusting the environment varible ` OMP_NUM_THREADS`.


The following pages describe the gory details to run FFEA simulations:

- @subpage ffea_iFile
- @subpage ffea_RB
- @subpage shortRange
- @subpage springs
- @subpage kineticApproach
- @subpage fmApproach
- @subpage ioFiles
- @subpage keywordReference

