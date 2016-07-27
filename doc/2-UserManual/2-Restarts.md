
Restarts and backups {#ffea_RB}
===============================




Restarts {#ffea_restarts}
=========================
The FFEA runner has the possiblity to continue or extend a simulation that 
 was completed or even crashed. When a simulation is to be continued,
 one will need:
  -  a trajectory, to get the positions of the nodes and, depending on
     the configuration, their velocities and forces. The new snapshots will be 
     appended after the last one found. 
  - a checkpoint file, that will indicate the state of the random number 
     generators, to ensure that correlations and consequent artifacts 
     never arise accidentally. 
  - a measurement file, that will be extended with the new measurements.


The previous run will have provided this set of files, their names being 
 specified in the FFEA input file through keywords:
  - ` trajectory_out_fname ` 
  - ` checkpoint_out ` 
  - ` measurement_out ` 
or having default values:
  - ` <ffea-input-file>.ftj ` 
  - ` <ffea-input-file>.fcp ` 
  - ` <ffea-input-file>.fm `  


In addition one needs to modify the already used FFEA input file, to indicate
 that the will to extend to simulations. Specifically one needs to set:
  - ` < restart = 1 >
  - ` < checkpoint_in = <valid-checkpoint-file> > ` 
while keeping ` trajectory_out_fname ` and ` measurement_out_fname ` unchanged.


There is a last rule that applies. Values for files ` checkpoint_out `
 and ` checkpoint_in ` must differ. Therefore, if in the previous run 
 ` checkpoint_out ` was not specified we recommed to 
 rename the old ` checkpoint_out ` ` .fcp ` file with its default name 
 to something else and use it as ` checkpoint_in ` value. 



Backups {#ffea_backups}
=======================
Output files will never be overwritten by the FFEA runner, but will be 
  renamed instead. This way, if a new simulation were launched by error 
  in the same folder where a previous simulation was run, every one 
  of the following files:
    * ` checkpoint_out ` 
    * ` trajectory_out_fname ` 
    * ` measurement_out_fname ` 

will be renamed to ` __<old-file-name>__bckp.<N> `, where ` N ` is the 
 smallest integer so that the resulting text does not match the name
 of any existing file in the working folder.



