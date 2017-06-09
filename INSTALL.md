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
       recent enough compiler. FFEA will compile with GCC 4.8 or Intel 15
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



