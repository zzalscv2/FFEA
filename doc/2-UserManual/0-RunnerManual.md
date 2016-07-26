FFEA runner {#userManual}
=========================

The FFEA runner is the program that will compute a FFEA trajectory given 
 an initial system under a set of conditions, all contained in the
 ` FFEA Input File `. The program is meant to run from the command line 
 typing:
 
      ffea <myInputFile.ffea>

It can run in parallel using multiple threads through OpenMP. By default, 
 it will try to use as many threads as cores found. One can control the 
 number this benaviour adjusting the environment varible ` OMP_NUM_THREADS`.


The following pages describe the gory details to run FFEA simulations:

- @subpage ffea_iFile
- @subpage ffea_RB
- @subpage shortRange
- @subpage kineticApproach
- @subpage fmApproach
- @subpage ioFiles
- @subpage keywordReference

