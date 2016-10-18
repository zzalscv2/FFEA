Running a Simulation {#ffeasimulationtut}
=============================

Before starting the simulation, let's have a look inside the script first.  
It's plain text, so open in your favourite text editor. You should see something like this:

![An example of an FFEA script file (.ffea)](ffeascript.png "GroEL FFEA Script")

The structure of file script tries to be as clear as possible, and should give you an idea of what is being calculated.
 Full detail on the algorithms and options available can be found in the [FFEA runner](\ref userManual) section,
  including the [input file reference](\ref keywordReference), which will also allow
   setting up more sophisticated simulations.
For demonstrating purposes, we need to change a few global parameters. 
 So, for now, go ahead and set `dt = 1e-15`, `num_steps = 1000` and `check = 10`. 
 This will set us up for a run of 1 nanosecond, output a frame every 10 simulation timesteps which are 0.01 nanoseconds in length.

Before starting the FFEA runner, you may want to consider setting your OMP_NUM_THREADS environment variable. The recommended number of threads is between 4 and 16. 16 will produce the best results, but will only be around 30-40% faster than running on 4 threads. For now, we will set

	export OMP_NUM_THREADS=4

Now, let's run ffea!

	ffea emd_5043_8ang.ffea

Once the simulation has finished (which should take less than 5 minutes), we can have a look at the structure in the [viewer](\ref FFEAviewertut), and also [perform some analysis](\ref FFEAanalysistut).
A small note; if a simulation crashes with an error talking about elements inverting, not to worry! The simulation timestep has been set just a little high. Set ` dt ` to a smaller value, and rerun.


