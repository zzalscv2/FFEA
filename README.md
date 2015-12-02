#
Install 
=======

This file describes how to install FFEA. The instructions in this file
are for the most common use cases, and cover the command line tools.


Prerequisites
=============

To install FFEA you need `cmake` and some third-party libraries:

   * CMake (>=2.8.11), the build system for FFEA.
     Required for building FFEA.
     https://cmake.org/

   * Boost (>=1.54.0). 
     Required compiled Boost library: program_options. 
     http://www.boost.org/

   * Eigen (>=3.2.2).
     FFEA uses Eigen within the Kinetics module.
     http://eigen.tuxfamily.org


Configuration
=============

To configure FFEA type:

    cmake 

within the FFEA source tree. It is generally advisable to build FFEA anywhere else. 
Therefore, we would recommend to:

    mkdir $FFEA_BUILD
    cd $FFEA_BUILD
    cmake $FFEA_SRC

where ` $FFEA_BUILD ` denotes the directory with the FFEA sources while 
  ` $FFEA_BUILD` is an arbitrary folder where the generated files will be placed.
Several options can be added to `cmake`, being the most important:

  * `-DCMAKE_INSTALL_PREFIX=<dir>`       -  (default /usr/local) installation directory
  * `-DCMAKE_BUILD_TYPE=<Debug|Release>` -  (default Release) build type
  * `-DCMAKE_CXX_COMPILER=<program>`     -  (default g++)  C++ compiler.

CMake will look for the required Boost and Eigen libraries. In the case they are not 
 installed in a standard place, you can help CMake either through 

  * configuring with ` -DCMAKE_PREFIX_PATH="Path-to-Eigen;Path-to-Boost" `,
  * exporting enviroment variables ` EIGEN_HOME `  and ` BOOST_ROOT ` to the corresponding 
      software folders
  * or configuring with ` -DEIGEN_HOME="Path-to-Eigen" ` and  ` -DBOOST_ROOT="Path-to-Boost `

Specific FFEA flags include:
  * `USE_FAST`    (default OFF) will try to find the best compiler flags in terms of performance.
  * `USE_OPENMP`  (default ON) will enable OpenMP parallel calculations.
  * `USE_OMP_MODE` (default 1) where:
    - 0 is for ` USE_OPENMP=OFF `.
    - 1 uses all the threads within blobs.
    - 2 uses one thread per blob.


Building
========
After configuring you will be able to build FFEA typing:

    make 

Optionally, you can build the documentation typing:

    make doc 

Finally, you can install FFEA:

    make install