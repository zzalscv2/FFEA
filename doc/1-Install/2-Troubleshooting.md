Troubleshooting {#troubleshooting}
=======
This tutorial will assume that the user does not have root access and must install packages locally. It's also assumed that the user will compile all of those packages from source. In the examples that follow, the software will be installed into $HOME/Software/LocalInstall, but any install folder can be used.
Before using this guide, you should already be familiar with the [basics of the Linux terminal](https://www.cheatography.com/davechild/cheat-sheets/linux-command-line/). You should also understand the [principles behind compiling, linking and building executables](https://www3.ntu.edu.sg/home/ehchua/programming/cpp/gcc_make.html)

## Installing the FFEA runner dependencies

### Download and install cmake
FFEA depends on cmake for building and testing. Download the latest cmake from cmake.org and extract it into a folder of your choosing. Open a terminal in that folder, or navigate to that folder in the terminal. There should be a file called 'bootstrap'. To execute it, type
```sh
./configure --prefix=$HOME/Software/LocalInstall
```
The --prefix flag configures the package to be installed into a non-standard location.
to build cmake, just type:
```sh
make
```
A cmake executable will be created. To install it, type
```sh
make install
```
Now that our software is in a non-standard location, we need to be able to tell other packages on the system where to look for it. To do this, we have to append the directory containing the cmake executable to our system path variable. If cmake is installed in the location above, we run:
```sh
export PATH=$PATH:$HOME/Software/LocalInstall/bin
```
in a Bash shell. These environment variables will reset every time you open a new terminal window. To preserve them, we can add them to our .bashrc file, which should be located in the home directory (if it's not, create a file called '.bashrc'). If you still can't see the file, remember to uncheck 'show hidden files' in your file manager, as by default, any file name beginning with a . is hidden.

The .bashrc file will execute all the commands listed in it when a new terminal is opened. For this reason, adding a line to the .bashrc file won't have any effect until the terminal is restarted.
### Download and install Boost
Boost is a general purpose C++ library with many useful features. Download the boost library from boost.org and extract the contents. Open a terminal window in the Boost folder. Boost comes with a bash script called bootstrap.sh. After extracting Boost, run
```sh
./bootstrap.sh --prefix=$HOME/Software/LocalInstall
```
Again, the --prefix flag tells the bootstrapper script that the package will be installed into a non-standard location.

When the bootstrapper has finished running, it will prompt you to run the script
```sh
./b2
```
to build Boost.

Finally, run
```sh
./b2 install
```
to install it.
If an error appears saying that there is no file called 'bzlib.h', you need to install the bzip2 library, which can be found at bzip.org. Download the source, extract it, navigate to that folder in the terminal, and run
```sh
make
```
and then
```sh
make install PREFIX=$HOME/Software/LocalInstall
```
to install it. Then, try installing Boost again from the very beginning. 
### Download and install Eigen
Eigen is a C++ template library for linear algebra.
Grab Eigen from eigen.tuxfamily.org.

Eigen's website and documentation state that Eigen doesn't need to be compiled, but for FFEA, it does. The easiest way to do so is to make a new folder in the Eigen folder called build_dir, switch to that folder, and inside it, run
```sh
cmake /path/to/Eigen/source -DCMAKE_INSTALL_PREFIX=$HOME/Software/LocalInstall
make install
```
You should then add Eigen's install folder to your CPATH environment variable. This will allow C and C++ compilers to find the Eigen header files.
If you are using a Bash shell, run this command:
```sh
export CPATH="~/Software/LocalInstall/usr/local/include/eigen3"
```
Finally, the ` EIGEN3_HOME ` environment variable will ensure that cmake can detect Eigen.
```sh
export EIGEN3_HOME="~/Software/LocalInstall/usr/local/include/eigen3"
```
Adding these commands to the .bashrc file is strongly recommended during the installation process of FFEA.

### Download and install Doxygen
Unless you are an FFEA developer, you can skip this step, as there is probably no reason to build the docs manually when they can be viewed at ffea.bitbucket.com.

Still, compiling and installing Doxygen in our local folder is as easy as:
```sh
./configure --prefix=$HOME/Software/LocalInstall
make
make install
```
However, Doxygen depends on Flex, which depends on GNU Bison. Building Flex is a bit of a pain, because Bison has its own environment variable that you'll need to change, or you may get an error saying that make can't find m4sugar.m4. Assuming you've installed Bison in the same fashion as previous packages, in this tutorial, modify the BISON_PKGDATADIR environment variable like so:
```sh
export BISON_PKGDATADIR=$HOME/Software/LocalInstall/share/bison
```
Flex and Doxygen should both now be able to build. Compiling the docs is one of the last steps, so for now, let's try and get FFEA itself to run.

## Installing FFEA
Finally, it's time to compile and install FFEA. Download and extract the FFEA source code. Make a new folder, separate from the source code, where you want the build files to live, and execute these commands:
```sh
cmake /path/to/ffea/source -DCMAKE_INSTALL_PREFIX=$HOME/Software/LocalInstall
make
make install
```
If you get an error saying that Boost or Eigen cannot be found, you may need to specify the locations of these manually. Although they can be specified by environment variables, the most reliable method is to use the following flags:
```sh
cmake /path/to/ffea/source -DEIGEN3_HOME=$HOME/Software/LocalInstall -DBOOST_ROOT=$HOME/Software/LocalInstall
make
```
where ` EIGEN3_HOME ` and ` BOOST_ROOT ` point to the folders where Eigen3 and Boost 
 were installed respectively.

If you encounter an internal compiler error, then you may try to get 
 around it by reducing the number of compiler optimisations. This will be achieved by 
 adding the -DUSE_FAST=OFF flag and running cmake again:
```sh
cmake /path/to/ffea/source -DUSE_FAST=OFF
```
USE_FAST tries to find the best compiler flags for performance, however it could fail
 for some platforms.

If you installed Doxygen, you can build the documentation from the same folder, using
```sh
make doc
```
And you can start the unit tests by running
```sh
make test
```
Finally, 
```sh
make install
```
will place the software in the desired location.


Now that you have installed FFEA, and assuming that the install folder is already 
 in your path, you should be able to run it just by typing ` ffea ` into the terminal. 
 If you get an error that says 
 ` error while loading shared libraries: libboost ... no such file or directory `, 
 then add the following line to your .bashrc file:
```sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Software/LocalInstall/lib
```
and restart the terminal.


## ffeatools
The ` ffeatools ` are a number of Python modules that are provided as part of a Python package. While 
 installing FFEA following the previous [notes](\ref installing-ffea) should have installed this package 
 under ` $HOME/Software/LocalInstall/lib/pythonX.Y/FFEA_python_modules `, this can be installed 
 as a standalone package, without the FFEA runner. In order to do that,
 open the FFEA source folder (containing ``setup.py``) in the terminal, and type

     python setup.py install

This will install the FFEA tools into your Python site-packages folder. Alternatively, a path
 can be specified through ` --prefix= `. 

FFEA tools use
 [NumPy](http://www.numpy.org/), so if you do not have root access onto your machine,
 or you do not have a package manager providing it, you may want to consider installing the [Anaconda Python distribution](https://www.continuum.io/downloads), which comes with a set of common scientific Python packages, and allows for packages to be installed and removed easily through a package manager (conda). Anaconda comes with a Python IDE called Spyder, which can be launched using the terminal

     spyder

In order to install, visit Anaconda's website to get a shell script, 
 which can be executed using the following command
```sh
bash Anaconda2-4.1.1-Linux-x86_64.sh
```
After some red tape, this will automatically install Anaconda.

