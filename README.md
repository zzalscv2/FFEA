Install 
=======

This file describes how to install FFEA. The instructions in this file
are for the most common use cases, and cover the command line tools.


Prerequisites
=============

To install FFEA you need:

   * C and a C++ compilers.   
     There is some C++ code written using 
       the C++11 standard, so CMake will ensure that you have a 
       recent enough compiler. Still, GCC 4.4 and Intel 13 have shown to work well. 

   * CMake (>=2.8.11).   
     Required for building FFEA.
     https://cmake.org/

and some third-party libraries:

   * Boost (>=1.54.0).   
     Required compiled Boost library: program_options, filesystem and system.
     http://www.boost.org/

   * Eigen (>=3.2.1).   
     FFEA uses Eigen within the Kinetics module.
     http://eigen.tuxfamily.org

      - Warning - All of the Eigen releases to date (<= 3.3~beta1-2) prove to 
          be incompatible with GCC >= 5 and the C++11 standard. Thus, to compile FFEA,
          you need to stick to either GCC 4, or use the latest development commit 
          from Eigen (3ba633e). 

 


   * Doxygen (>= 1.8) [OPTIONAL]   
     It will be used to build the documentation. http://www.doxygen.org



Configuration
=============

It is generally advisable to configure and compile FFEA outside of the source tree. 
Therefore, to configure FFEA, we would recommend to:

    mkdir $FFEA_BUILD
    cd $FFEA_BUILD
    cmake $FFEA_SRC

where ` $FFEA_SRC ` denotes the directory with the FFEA sources while 
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

Optionally, if Doxygen was found at configure time, 
 you can build the documentation typing:

    make doc 

there are some mathematical formulae that will not render correctly if latex and ghostview are not found. You can read the documentation with a browser, if firefox was the browser available to you the coomand would be:

firefox $FFEA_HOME/share/doc/html/index.html &


Finally, you can install FFEA:

    make install

You can run the software with the command:

ffea <input-file.ffea>


The only caveat in this last step is that <input-file.ffea> points to a number of input files. Details of all these files are in the documentation. You must be in the right folder so that the relative paths in <input-file.ffea> make sense.