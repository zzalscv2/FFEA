FFEA runner {#userManual}
=========================

The FFEA runner is the program that will compute a FFEA trajectory given 
 an initial system under a set of conditions, as described in the
 ` FFEA Input File `. The program is compiled in two flavours: ` ffea ` 
 and ` ffea_mb ` and is meant to run from the command line typing:
 
      ffea <myInputFile.ffea>

Both ` ffea ` and ` ffea_mb ` are parallel programs using OpenMP. The 
 former, ` ffea ` uses all the threads within every blob, while the later, 
 ` ffea_mb ` assigns one (or more) blobs to every thread. Ideally, 
 ` ffea ` is to be used in systems where there is a single enormous blob
 (in terms of number of nodes), while ` ffea_mb ` performs much better in 
 a system with many blobs.


A number of optional arguments are available:

	'-h/--help' - Prints out a help text to the command prompt
	'-v/--verbose' - Prints out to the command line a user defined level of information
	'-d/--no-detail' - Only writes measurement data to files for the total global system, not each individual blob and pair of blobs (the .fdm file is not produced)
	'-m/--mode' - Select the mode of FFEA you would like to run: Full simulation [default] (0), Linear Elastic Model (1), Linear Dynamic Model (2) or Timestep Calculator (3). 
 
FFEA can run in parallel using multiple threads through OpenMP. By default, 
 it will try to use as many threads as cores found. One can control the 
 number this behaviour adjusting the environment variable ` OMP_NUM_THREADS`.


The following pages describe the gory details to run FFEA simulations:

- @subpage ffea_iFile
- @subpage ffea_RB
- @subpage shortRange
- @subpage restraints
- @subpage kineticApproach
- @subpage fmApproach
- @subpage ctforces
- @subpage ioFiles
- @subpage keywordReference

