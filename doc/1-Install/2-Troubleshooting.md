Troubleshooting {#troubleshooting}
=======
Before using this guide, you should already be familiar with the [basics of the linux terminal](https://www.cheatography.com/davechild/cheat-sheets/linux-command-line/). You should also understand the [principles behind compiling, linking and building executables](https://www3.ntu.edu.sg/home/ehchua/programming/cpp/gcc_make.html)

This tutorial will assume that the user does not have root access and must install packages locally. It's also assumed that the user will compile all of those packages from source. In the examples that follow, the software will be installed into $HOME/Software/LocalInstall, but any install folder can be used.

## Installing the FFEA runner

### Download and install cmake
FFEA depends on cmake for building and testing. Download the latest cmake from cmake.org and extract it into a folder of your choosing. Open a terminal in that folder, or navigate to that folder in the terminal. There should be a file called 'bootstrap'. To execute it, type
```sh
./bootstrap
```
to build cmake, just type:
```sh
make
```
A cmake executable will be created. To install it, type
```sh
make install
```
If you don't have root access on your machine, you won't be able to install like this, as you don't have access to the folders it wants to install into. The following command will change the installation directory, so that cmake installs into the home folder:
```sh
export DESTDIR="$HOME/Software/LocalInstall"
```
Here we have used the folder /Software/LocalInstall, but it can go anywhere.

Now that our software is in a non-standard location, we need to be able to tell other packages on the system where to look for it. To do this, we have to append the directory containing the cmake executable to our system path variable. If cmake is installed in the location above, we run:
```sh
export PATH=$PATH:$HOME/Software/LocalInstall/usr/local/bin
```
These environment variables will reset every time you open a new terminal window. To preserve them, we can add them to our .bashrc file, which should be located in the home directory (if it's not, create a file called '.bashrc'). If you still can't see the file, remember to uncheck 'show hidden files' in your file manager, as by default, any filename beginning with a . is hidden.

The .bashrc file will execute all the commands listed in it when a terminal is opened. For this reason, adding a line to the .bashrc file won't have any effect until the terminal is restarted.
### Download and install Boost
Boost is a general purpose C++ library with many useful features. Download the boost library from boost.org and extract the contents. Open a terminal window in the Boost folder. Boost comes with a bash script called bootstrap.sh. After extracting Boost, run
```sh
./bootstrap.sh --prefix=$HOME/Software/LocalInstall
```
The --prefix flag tells the bootstrapper script that the package will be installed into a non-standard location.

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
to install it. Then, try installing Boost again. Once you're done instaling Boost, don't delete the source code. Programs built with Boost depend on the source code, headers and libraries, which are all in different locations.
### Download and install Eigen
Eigen is a C++ template library for linear algebra.
Grab Eigen from eigen.tuxfamily.org.

Eigen's website and documentation state that Eigen doesn't need to be compiled, but for FFEA, it does. The easiest way to do so is to make a new folder in the Eigen folder called build_dir, switch to that folder, and inside it, run
```sh
cmake /path/to/seigen/source
make install
```
You should then add Eigen's install folder to your CPATH environment variable. This will allow C and C++ compilers to find the Eigen header files.
If you are using a Red Hat based distro (e.g. CentOS, Red Hat, Fedora, or Scientific Linux), run this command:
```sh
export CPATH="~/Software/LocalInstall/usr/local/include/eigen3"
```
The directory structure of Debian-based distros (e.g. Ubuntu and Debian) is different, so you should only run:
```sh
export CPATH="~/Software/LocalInstall/usr/local/include"
```
Finally, the EIGEN_HOME environment variable will ensure that cmake can detect Eigen.
```sh
EIGEN_HOME="~/Software/LocalInstall/usr/local/include/eigen3"
```
Adding these commands to the .bashrc file is strongly recommended.

If cmake still won't recognize eigen, try entering the cmake folder from the eigen source, finding the file named 'FindEigen3.cmake' and copying it to the cmake folder in the FFEA folder.
### Download and install Doxygen
Unless you are an FFEA developer, you can skip this step, as there is probably no reason to build the docs manually when they can be viewed at ffea.bitbucket.com.

Doxygen depends on Flex, which depends on GNU Bison. Building Flex is a bit of a pain, because Bison has its own environment variable that you'll need to change, or you may get an error saying that make can't find m4sugar.m4. Assuming you've installed Bison in the same fashion as previous packages, in this tutorial, modify the BISON_PKGDATADIR environment variable like so:
```sh
export BISON_PKGDATADIR=$HOME/Software/LocalInstall/usr/local/share/bison
```
Flex and Doxygen should both now be able to build. Compiling the docs is one of the last steps, so for now, let's try and get FFEA itself to run.
### Installing FFEA
Finally, it's time to compile and install FFEA. Download and extract the FFEA source code. Make a new folder, separate from the source code, where you want the build files to live, and execute these commands:
```sh
cmake /path/to/ffea/source
make
make install
```
If you get an error saying that Boost or Eigen cannot be found, you may need to specify the locations of these manually. Although they can be specified by environment variables, the most reliable method is to use the following flags:
```sh
cmake /path/to/ffea/source -DEIGEN_HOME=$HOME/Software/src/eigen-eigen-dc6cfdf9bcec -DBOOST_ROOT=/localhome/username/Software/src/boost_1_61_0 -DBOOST_INCLUDEDIR=/localhome/username/Software/LocalInstall/include/boost -DBOOST_LIBRARYDIR=/localhome/username/Software/LocalInstall/lib -DCMAKE_MODULE_PATH=/localhome/py12rw/Software/src/eigen-eigen-dc6cfdf9bcec/cmake
make
make install
```
BOOST_LIBRARYDIR and BOOST_INCLUDEDIR are two directories created when installing Boost, and will be wherever you specified Boost to be installed, include in /include/boost and library in /lib. BOOST_ROOT is wherever you left the source code for Boost. Similarly, EIGEN_HOME is the source for Eigen. Finally, CMAKE_MODULE_PATH is the path to eigen's cmake modules that we discussed earlier.

If you encounter an internal compiler error, then you may be using an old version of gcc (for example, 4.4 has this issue). If you cannot update gcc, you may be able to get around it by reducing the number of compiler optimisations. There are a few flags that have this function:
```sh
cmake /path/to/ffea/source -DUSE_FAST=OFF -DUSE_OPENMP=OFF
```

USE_FAST tries to find the best compiler flags for performance. USE_OPENMP enables openMP parallel calculations. Try adding the -DUSE_FAST=OFF flag and running cmake again.

Now that you have installed FFEA, you should be able to run it just by typing 'ffea' into the terminal. If you get an error that says 'error while loading shared libraties: libboost ... no such file or directory', then add the following line to your .bashrc file:
```sh
export LD_LIBRARY_PATH="/localhome/yourusername/Software/LocalInstall/lib"
```
and restart the terminal.

If you installed doxygen, you can build the documentation from the same folder, using
```sh
make doc
```
And you can start the unit tests by running
```sh
make test
```
## ffeatools
ffeatools have been tested on Python versions after 2.6x and 2.7x. If you do not already have Python, the easiest way to install it is probably through Anaconda, available at continuum.io/downloads. The Anaconda distribution will install a selection of common Python packages, a package manager (conda) and a basic Python IDE (spyder). Even if you do already have Python, it is recommended, and it will not overwrite the default python (although it will launch as 'Python' from the terminal).

Anaconda's website provides a shell script, which can be executed using the following command
```sh
bash Anaconda2-4.1.1-Linux-x86_64.sh
```
After some red tape, this will automatically install Anaconda.

The FFEA tools are installed alongisde the runner, and can be accessed using the following command:
```sh
ffeatools
```
