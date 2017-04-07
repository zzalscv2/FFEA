Overview {#overview}
=========

Fluctuating Finite Element Analysis is a new molecular modelling technique, built from the ground-up to support systems that are larger and more complex than those modelled by atomistic molecular dynamics. Instead of modelling biological systems as a collection of connected atoms, it models them as 3D volumes comprised of tetrahedrons. Unlike previous coarse-grained models, the models FFEA generates are visco-elastic continuum solids. Unlike other applications of Finite Element Analysis, these systems are subject to thermal fluctuations.

This technique has the potential to model large, complex systems, made of many molecules, and complex processes at the frontiers of molecular biology. As it does not not require an atomistic level of detail, it can also be used to simulate biological molecules that cannot be imaged using X-ray crystallography.


Features  {#features}
========

 * Protein Interactions:
  * A 6-12 Lennard-Jones Potential.
  * A repulsive potential that is proportional 
        to the overlapping volume.
  * Specific interactions defined using precomputed potentials.
        More documentation can be found [here](\ref fmApproach).
  * Coulombic interactions [EXPERIMENTAL].
 * Kinetic state changes can be simulated together with the continuum model to
    account for conformational changes and binding events.
 * Conversion tools for EM density data and atomistic structures into FFEA simulations.
 * A plugin for PyMOL, allowing the visualisation of FFEA systems and trajectories.
 * Analysis tools (equilibration, Euler characteristic, principal component analysis, geometric measurements) available on the command line and under a Python API.
 * Extensive test suite including checks of FFEA's simulation output against analytical results.



Publications  {#publications}
============

   * Methodology
       * Oliver R., Read D. J., Harlen O. G. & Harris S. A. ["A Stochastic finite element model for the dynamics of globular macromolecules"](http://www.sciencedirect.com/science/article/pii/S0021999112007589) (2013) J. Comp. Phys. 239, 147-165.
       * Patargias G. N., Harris S. A. & Harding J. ["A demonstration of the inhomogeneity of the local dielectric response of proteins by molecular dynamics simulations."](https://www.ncbi.nlm.nih.gov/pubmed/20572740) (2010) J. Chem. Phys. 132, 235103.
   * Applications
       * Richardson R., Papachristos K., Read D. J., Harlen O. G., Harrison M. A., Paci E., Muench S. P. & Harris S. A ["Understanding the apparent stator-rotor connections in the rotary ATPase family using coarse-grained computer modelling"](https://www.ncbi.nlm.nih.gov/pubmed/25174610) (2014), Proteins: Struct., Funct., Bioinf., 82, 3298-3311.
   * Reviews
       * Gray A., Harlen O. G., Harris S. A., Khalid S., Leung Y. M., Lonsdale R., Mulholland A. J., Pearson A. R., Read D. J. & Richardson R. A. ["In pursuit of an accurate spatial and temporal model of biomolecules at the atomistic level: a perspective on computer simulation"](https://www.ncbi.nlm.nih.gov/pubmed/25615870), Acta Cryst. (2015) D71, 162-172.
       * Oliver R. , Richardson R. A., Hanson B., Kendrick K., Read D. J., Harlen O. G. & Harris S. A. ["Modelling the Dynamic Architecture of Biomaterials Using Continuum Mechanics"](http://link.springer.com/chapter/10.1007%2F978-3-319-09976-7_8), Protein Modelling, G. Náray-Szabó, Editor. (2014) Springer International Publishing. p. 175-197.
       * Hanson B., Richardson R., Oliver R., Read D. J., Harlen O. & Harris S. ["Modelling biomacromolecular assemblies with continuum mechanics"](https://www.ncbi.nlm.nih.gov/pubmed/25849915) Biochem. Soc. Trans. (2015), 43, 186-192.


Technology  {#technology}
============
 
   * [Boost](http://www.boost.org) (>=1.54.0) is used 
     for ease of programming 
     at the initialisation phase. Modules "system", "filesystem" and 
     "program-options" are required. A bundle of version 1.63 is shipped with FFEA.
   * [Eigen](http://eigen.tuxfamily.org) (>=3.2.1).   
     FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models. CMake will download and use Eigen 3.3.2 if not told otherwise.
   * [RngStreams](http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/)<sup>[1](#RngStreams1)</sup><sup>,[2](#RngStreams2)</sup>
        is shipped with FFEA and used as Random Number Generator (RNG). RngStreams 
        allows the FFEA to **safely** generate random numbers when running 
        on a number of threads, as well as safe restarts, recovering the state 
        of the RNGs in the last saved time step.
   * [Tet_a_tet](https://github.com/erich666/jgt-code/blob/master/Volume_07/Number_2/Ganovelli2002/tet_a_tet.h)<sup>[3](#tetatetpaper)</sup>
        is shipped with FFEA and used to detect element overlapping 
        in the steric repulsion module. 
   * [Doxygen](http://www.doxygen.org) (>= 1.8) [OPTIONAL] 
        is used to generate the documentation. 
   * [PyMOL](https://www.pymol.org) (>=1.8) [OPTIONAL] can 
        be used, using the plugin we provide,
        to visualise FFEA systems and trajectories
        as well as molecular and EM systems.
   * [mtTkinter](http://tkinter.unpythonic.net/wiki/mtTkinter) (0.4) is shipped 
        with FFEA and used in the PyMOL plugin, allowing safe threading. 
   * [GTS](http://gts.sourceforge.net) (>=0.7.6)[OPTIONAL]. The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.
   * [NETGEN](https://sourceforge.net/projects/netgen-mesher/) 
   or [TETGEN](http://wias-berlin.de/software/tetgen/) [OPTIONAL]. 
     Programs which convert surface profile into volumetric meshes 
        to be used by FFEA.
   * [pyPcazip](https://pypi.python.org/pypi/pyPcazip) [OPTIONAL]
     Some of the Python FFEA analysis tools interact with these 
     Principal Component Analysis library in order to generate the standard
     PCA output (eigensystems, projections, animations etc)
     obtained from standard from equivalent MD simulations.


<a name="RngStreams1">1</a>: P L'Ecuyer, "Good Parameter Sets for Combined Multiple Recursive Random Number Generators", (1999) Oper Res, 47(1):159-164. <br> 
<a name="RngStreams2">2</a>: P L'Ecuyer et al. "An Objected-Oriented Random-Number Package with Many Long Streams and Substreams", (2002) Oper Res, 50(6):1073-1075. <br>
<a name="tetatetpaper">3</a>:  F Ganovelli, et al. "Fast tetrahedron-tetrahedron overlap algorithm", (2002) J.Graphics Tools. 7(2):17-25.



Contribute
==========

FFEA is maintained by a small but dedicated team at the University of Leeds. If you want to see where we can take FFEA, then you can:

   * Use the software for something cool
   * Send bug reports and feature requests to our [issue tracker](https://bitbucket.org/sohpc-ffea/ffea/issues)
   * [Fork us](https://bitbucket.org/sohpc-ffea/ffea/fork)


FFEA Team  {#FFEAteam}
==========

### Code: ###
   * Albert Solernou
   * Ben Hanson
   * Robin Richardson
   * [Rob Welch](http://robwel.ch/)


### Theory: ###
   * [Oliver Harlen](https://www.maths.leeds.ac.uk/index.php?id=263&uid=1025)
   * [Sarah Harris](http://www.comp-bio.physics.leeds.ac.uk/)
   * [Daniel Read](http://www1.maths.leeds.ac.uk/~djread/)
   * Robin Oliver
   * Robin Richardson
   * Ben Hanson
   * Albert Solernou


Thanks {#thanks}
=======
We want to thank everybody who has helped in making FFEA possible, from 
 summer students, to experimental professors, we would not be here without you:
   * Stan Burgess
   * Stephen Muench
   * Kerrie Smith
   * Joanna Leng
   * Thijs van der Heijden
   * Kees Storm
   * Paul van der Schoot 
   * Toni Collis
   * Neelofer Banglawala
   * Jana Boltersdorf
   * Ondřej Vysocký
   * Guanhao Lu
   * Jonathan Boyle
   * Mike Croucher
   * Christopher Woods


Installation {#install}
============

This document gives instructions to build and install the FFEA package,
 consisting of the FFEA runner and the FFEA tools. 


Prerequisites {#prerequisites}
=============

To install FFEA you need:

   * C and C++ compilers.   
     There is some C++ code written using 
       the C++11 standard, and so CMake will ensure that you have a 
       recent enough compiler. FFEA will compile with GCC 4.4 or Intel 15
       but we recommend to use an up to date compiler to get the best performance.

   * [CMake](https://cmake.org) (>=2.8.11).   
     Required for building FFEA.

   * [Doxygen](http://www.doxygen.org) (>= 1.8) [OPTIONAL]   
     It will be used to build the documentation.

   * [Python](https://www.python.org/) (>= 2.6), is used to run some tests to verify that the FFEA runner was correctly built, as well as in the FFEA tools.

   * [NumPy](http://www.numpy.org/) is needed to run some tests as well as in 
        several FFEAtools.


FFEA uses Boost and Eigen. To make your life easier, **the code is shipped whih 
 a subset of Boost (v. 1.63), and Eigen (v 3.3.2)** will be downloaded
 at configure time. Still, you are welcome to use your own versions of the libraries.

   * [Boost](http://www.boost.org) (>=1.54.0) [OPTIONAL]   
     FFEA uses modules: program_options, filesystem and system.

   * [Eigen](http://eigen.tuxfamily.org) (>=3.2.1) [OPTIONAL]
   FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models.

   > Warning - GCC >= 5 will require version 3.2.10 or higher for Eigen. Earlier versions (including 3.3 release candidates, marked internally as 3.2.91 and higher) did prove to be incompatible with GCC >= 5 and using C++11 standard. You may also need a newer version of Boost. 
     

Configure {#configure}
=========

It is generally advisable to configure and compile FFEA outside of the source tree. 
Therefore, to configure FFEA, we would recommend to:

    mkdir $FFEA_BUILD
    cd $FFEA_BUILD 
    cmake $FFEA_SRC [OPTIONS]

where ` $FFEA_SRC ` denotes the directory with the FFEA sources while 
  ` $FFEA_BUILD` is an arbitrary folder where the generated files will be placed.
Several ` [OPTIONS] ` can be added to `cmake`, being the most important ones:

  * `-DCMAKE_INSTALL_PREFIX=<install_dir>`       -  (default /usr/local) installation directory
  * `-DCMAKE_BUILD_TYPE=<Debug|Release>` -  (default Release) build type
  * `-DCMAKE_CXX_COMPILER=<program>`     -  (default g++)  C++ compiler.

By default, CMake will use the subset of Boost we bundled, and will download 
 Eigen 3.3.2. If you want to use your own installations you can still do so
 through:
  
  * `-DUSE_BOOST_INTERNAL=<ON|OFF>` - default ` ON `
  * `-DUSE_EIGEN3_INTERNAL=<ON|OFF>` - default ` ON `
 
If you decide to do so, CMake will look for the required Boost and/or Eigen libraries.
 In the case they are not 
 installed in a standard place, you can help CMake either through: 

  * ` -DCMAKE_PREFIX_PATH="Path-to-Eigen;Path-to-Boost" `,
  * ` -DEIGEN_HOME="Path-to-Eigen" ` and  ` -DBOOST_ROOT="Path-to-Boost" `
  * or exporting environment variables ` EIGEN_HOME `  and ` BOOST_ROOT ` 
      to the corresponding software folders

Additional specific FFEA flags include:

  * `USE_FAST`    (default ON) will try to find the best compiler flags in terms of performance.
  * `USE_OPENMP`  (default ON) will enable OpenMP parallel calculations.
  * `USE_OMP_MODE` (default 1) where:

    - 0 is for ` USE_OPENMP=OFF `.
    - 1 uses all the threads within blobs. The only option to get OMP parallism if simulating a single body. 
    - 2 uses one thread per blob. Performs significantly better if there is more than a single body.

Thus, for production runs, one could configure the package typing:

    cmake $FFEA_SRC -DCMAKE_INSTALL_PREFIX=$HOME/softw/ffea






Build {#build} 
=====
After configuring you will be able to build FFEA typing:

     make 

Optionally, if Doxygen was found at configure time, 
 you can build the documentation typing:

     make doc 

There are some mathematical formulae that will not render correctly if latex
  and ghostview are not found.


Finally you may want to check your installation running a provided suite of tests, 
 either sequentially:
  
     make test

or concurrently:

     ctest -j <number-of-processes> 


  

Install {#install} 
=======
The following command will install FFEA:

    make install

either to the folder specified through ` -DCMAKE_INSTALL_PREFIX `
  or into a default folder (where you may need administrative privileges).


The FFEA_runner, ` ffea `, as well as the ffeatools, ` ffeatools ` will be found 
 in ` $FFEA_HOME/bin `. Instructions on how to use them can be read 
 [here](\ref userManual) and [here](../../ffeamodules/html/index.html) respectively. 
 In addition, the ` ffeatools ` Python package will be found in 
 ` $FFEA_HOME/lib/python<version>/site-packages `

If you built the documentation you will be able to read it with a web browser, 
  and so if firefox was the browser of your choice, and you installed 
  FFEA in ` $FFEA_HOME `, the command would be:

      firefox $FFEA_HOME/share/ffea/doc/html/index.html &



Finally, a plugin to visualise systems and trajectories in 
 [PyMOL](https://www.pymol.org) should be found in:

     $FFEA_HOME/share/ffea/plugins/pymol/ffea.tar.gz


In order to use it, one would need to run PyMOL (>=1.6), and then click on
  ` Plugin ` -> ` Plugin Manager `, and on the new window, go to tab 
  ` Install New Plugin `, click ` Choose file... ` and finally find and 
  select ` ffea.tar.gz ` from your disk. You will be asked to install the 
  plugin either into a local folder or a global folder. If the 
  folder does not exist, PyMOL will ask your permission on creating the folder, 
  and will install the plugin properly. However, just in this case,
  it but will bump ` Plugin FFEAplugin has
  been installed but initialization failed `, and you'll need to restart PyMOL
  to use the plugin.



Working environment {#workingEnvironment}
===================

Executing ` ffea ` and ` ffeatools ` is probably what most users will wish, so 
 UNIX users may find convenient to add the install folder in the ` PATH `:

      export PATH=$FFEA_HOME/bin:$PATH

In addition, in order to have direct access to the python modules that integrate 
 the FFEA tools, and specially to those users willing to write new measure tools,
 the ` PYTHONPATH ` environment variable should be also updated: 

     export PYTHONPATH=$FFEA_HOME/lib/python<version>/site-packages



Useful packages
===============

Once FFEA has been installed, users may want to provide themselves with some
 extra packages that have proved to be useful at setting up the system
 to simulate, as well as at analysing the results:

   * [PyMOL](https://www.pymol.org) (>=1.7) can
        be used to visualise FFEA systems and trajectories
        once the plugin has been installed.


   * [GTS](http://gts.sourceforge.net) (>=0.7.6)[OPTIONAL]. The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.


   * [NETGEN](https://sourceforge.net/projects/netgen-mesher/)
   or [TETGEN](http://wias-berlin.de/software/tetgen/) [OPTIONAL].
     Programs which convert surface profile into volumetric meshes
        to be used by FFEA.


   * [pyPcazip](https://pypi.python.org/pypi/pyPcazip) [OPTIONAL]
     Some of the Python FFEA analysis tools interact with these
     Principal Component Analysis library in order to generate the standard
     PCA output (eigensystems, projections, animations etc)
     obtained from standard from equivalent MD simulations.

Some notes on how to use these tools in relation to FFEA can be found
 in the [Tutorial](\ref Tutorial). However, mastering these tools
 may imply consulting the documentation provided by the packages themselves.



