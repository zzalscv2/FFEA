Quickstart {#quickstart}
============

Prerequisites
=============

To install FFEA you need:

   * C and C++ compilers.   
     There is some C++ code written using 
       the C++11 standard, so CMake will ensure that you have a 
       recent enough compiler. Still, GCC 4.4 and Intel 13 have shown to work well. 

   * [CMake](https://cmake.org) (>=2.8.11).   
     Required for building FFEA.

   * [Boost](http://www.boost.org) (>=1.54.0).   
     Required compiled Boost library: program_options, filesystem and system.

   * [Eigen](http://eigen.tuxfamily.org) (>=3.2.1). FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models.

   > Warning - version 3.2.10 will be needed if using GCC >= 5. Earlier versions did prove to be incompatible with GCC >= 5 and the C++11 standard. 
     
   * [Doxygen](http://www.doxygen.org) (>= 1.8) [OPTIONAL]   
     It will be used to build the documentation.

   * [Python](https://www.python.org/) (>= 2.6), is used to run some tests to verify that the FFEA runner was correctly built, as well as in the FFEA tools.

   * [NumPy](http://www.numpy.org/) is needed at run time in some FFEA tools.


Configure
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

CMake will look for the required Boost and Eigen libraries. In the case they are not 
 installed in a standard place, you can help CMake either through: 

  * configuring with ` -DCMAKE_PREFIX_PATH="Path-to-Eigen;Path-to-Boost" `,
  * exporting environment variables ` EIGEN_HOME `  and ` BOOST_ROOT ` to the corresponding 
      software folders
  * or configuring with ` -DEIGEN_HOME="Path-to-Eigen" ` and  ` -DBOOST_ROOT="Path-to-Boost" `

Additional specific FFEA flags include:

  * `USE_FAST`    (default ON) will try to find the best compiler flags in terms of performance.
  * `USE_OPENMP`  (default ON) will enable OpenMP parallel calculations.
  * `USE_OMP_MODE` (default 1) where:

    - 0 is for ` USE_OPENMP=OFF `.
    - 1 uses all the threads within blobs.
    - 2 uses one thread per blob.

Thus, for production runs, one could configure the package typing:

    cmake $FFEA_SRC -DCMAKE_INSTALL_PREFIX=$HOME/softw/ffea -DUSE_FAST=ON




Build
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


  

Install
=======
The following command will install FFEA:

    make install

either to the folder specified through ` -DCMAKE_INSTALL_PREFIX `
  or into a default folder.


If you built the documentation you will be able to read it wit a browser, 
  and so if firefox was the browser available to you, and you installed 
  FFEA in ` $FFEA_HOME `, the command would be:

      firefox $FFEA_HOME/share/ffea/doc/html/index.html &


The FFEA_runner, ` ffea `, as well as the ffeatools, ` ffeatools ` will be found 
 in ` $FFEA_HOME/bin `. Instructions on how to use them can be read 
 [here](\ref userManual) and [here](\ref FFEAtools) respectively. 

In addition, a plugin to visualise systems and trajectories in 
 [PyMOL](https://www.pymol.org) should be found in:

     $FFEA_HOME/share/ffea/plugins/pymol/ffea.tar.gz


In order to use it, one would need to run PyMOL (>=1.8), and then click on
  ` Plugin ` -> ` Plugin Manager `, and on the new window, go to tab 
  ` Install New Plugin `, click ` Choose file... ` and finally find and 
  select ` ffea.tar.gz ` from your disk.



Working environment {#workingEnvironment}
===================

Executing ` ffea ` and ` ffeatools ` is probably what most users will wish, so 
 UNIX users may find convenient to add the install folder in the ` PATH `:

      export PATH=$FFEA_HOME/bin:$PATH

In addition, in order to have access to the python modules that integrate 
 the FFEA tools, and specially to those users willing to write new measure tools,
 the ` PYTHONPATH ` environment variable should be also updated: 

     export PYTHONPATH=$FFEA_HOME/lib/python<version>/FFEA_python_modules

