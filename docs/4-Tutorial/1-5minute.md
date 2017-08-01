FFEA's Famous 5-Minute Tutorial {#fiveminute}
=========================

If you haven't already, [install FFEA!](\ref install) This shouldn't take too much of your alotted 5 minutes, all dependencies are automatically taken care of by cmake.

Once FFEA is installed, pick a protein that you want to simulate. Your protein should be globular (e.g. no long stringy bits - those are supported but require extra preparation) and large (ten nanometers long\wide at least!). Ideally, you should download your protein from [the EMDB](https://www.ebi.ac.uk/pdbe/emdb/), but you can also use the [regular PDB](http://www.rcsb.org/pdb/home/home.do). If you're not sure to pick, we like to [use GroEL](http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-5043) as an example.

Open up your molecule in a viewer (we recommend [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/)). Under the 'visual analysis' link on the EMDB, you should see a 'recommended contour level' setting. Set this contour level in Chimera, and check that the model looks right - you should see a smooth, consistent, connected structure. If not, adjust the contour level.

FFEA is split into two parts. The FFEA simulation itself is simply called 'ffea'. The other part is 'ffeatools', which contains the necessary tools to prepare, analyse and visualise simulations. The ffeatools package contains a tool called 'automodel', which will automatically perform all of the basic model generation features.

To run automodel, open a terminal in the same folder as your downloaded .map file, and run the command

```
ffeatools automodel [isolevel] [filename]
```

You can also run `ffeatools automodel -h` to see some of automodel's other parameters, like the model granularity or material parameters (density, stiffness, etc). They're not neccessary for testing FFEA, but setting them correctly is critical to producing a physically-accurate FFEA model.

After FFEA automodel has run, it will prompt you to run an FFEA simulation on the model you've just created, or edit the FFEA script file. If you want your trial simulation to run a bit faster, open up the .ffea script file in your favourite text editor, and set `dt = 1e-15`, `num_steps = 1000` and `check = 10`. If you're curious, [take a look at a more detailed breakdown of the FFEA input format](\ref ffea_iFile).

To run an ffea simulation, first enter the command
```
export OMP_NUM_THREADS=4
```
This will tell OpenMP to use four CPU threads. Then type
```
ffea [your .ffea script filename]
```
To run FFEA. FFEA should load up and start producing output immediately, and then terminate after 1000 steps and 100 frames.

That's it! You've now run your first simulation. To take a look at the results, install PyMOL (you could do so using [Anaconda](https://www.continuum.io/downloads) and `conda install pymol`), [install the FFEA viewer](\ref FFEAanalysistut), and open up the .ffea script file.

This tutorial, and FFEA automodel, only represent the most basic FFEA features. To see more, please [check out the full tutorial](\ref Tutorial), and the [FFEA runner documentation](\ref userManual).

