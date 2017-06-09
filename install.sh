#!/bin/bash

#wicked and unnessecary ascii art logo
echo ""
echo " ███████╗███████╗███████╗ █████╗"
echo " ██╔════╝██╔════╝██╔════╝██╔══██╗"
echo " █████╗  █████╗  █████╗  ███████║"
echo " ██╔══╝  ██╔══╝  ██╔══╝  ██╔══██║"
echo " ██║     ██║     ███████╗██║  ██║"
echo " ╚═╝     ╚═╝     ╚══════╝╚═╝  ╚═╝"
echo "FLUCTUATING FINITE ELEMENT ANALYSIS"
echo ""
              

# Check for mandatories              

if ! type "g++" > /dev/null && ! type "clang" > /dev/null && ! type "icc" > /dev/null; then
    echo "Error: no c++ compiler. Please install one!"
    exit 1
fi

if ! type "cmake" > /dev/null; then
    echo "Error: no cmake. Please install cmake before continuing!"
    exit 1
fi
                    
#Check for optional prerequisites
tetgen_exists=true
pymol_exists=true
doxygen_exists=true
numpy_exists=true

if ! type "tetgen" > /dev/null; then
    tetgen_exists=false
    echo "Warning: tetgen not found. FFEA can run without tetgen, but you will not be able to generate models."
fi

if ! type "doxygen" > /dev/null; then
    doxygen_exists=false
    echo "Warning: doxygen not found. FFEA will run, but you will not be able to compile the documentation."
fi

if ! type "pymol" > /dev/null; then
    pymol_exists=false
    echo "Warning: pymol not found. FFEA will run, but you will not be able to visualize the results."
fi

if ! python -c 'import numpy;'; then
    echo "Warning: numpy not found. Most of ffeatools  will not work."
    numpy_exists=false
fi

if ! $tetgen_exists || ! $doxygen_exists || ! $pymol_exists || ! $numpy_exists; then
    echo "Some optional dependencies were not found. Continue installing? (y\n):"
    read $missing_deps
    if [ "$missing_deps" == "n" ]; then
        exit 1
    fi
fi



#Ask user for path
echo "Where to install? Hit return without typing anything to install in the default directory (~/ffea)"
home=~
read install_dir

if [[ "$install_dir" == "" ]]; then
    install_dir=$home"/ffea"
fi

#Do out of source build

echo "Compiling and installing FFEA... (this might take a few minutes)"

mkdir ffea_src
shopt -s extglob # turn on extglob so that globs will work
mv !(ffea_src|install.sh) ffea_src # move everything except this script to the ffea_src folder
mkdir ffea_build
cd ffea_build
cmake ../ffea_src -DCMAKE_INSTALL_PREFIX=$install_dir || exit 1
make || exit 1
make install || exit 1

#Add to the .bashrc
path_export='export PATH="$PATH:'
path_export+=$install_dir
path_export+='/bin"'

ffeatools_path_export='export PATH="$PATH:'
ffeatools_path_export+=$install_dir
ffeatools_path_export+='/lib/python2.7/site-packages/ffeatools"'ls

ffeatools_pythonpath_export='export PYTHONPATH="$PYTHONPATH:'
ffeatools_pythonpath_export+=$install_dir
ffeatools_pythonpath_export+='/lib/python2.7/site-packages/ffeatools/modules"'

echo "#FFEA stuff" >> ~/.bashrc

echo $path_export >> ~/.bashrc
echo $ffeatools_path_export >> ~/.bashrc
echo $ffeatools_pythonpath_export >> ~/.bashrc

eval $path_export
eval $ffeatools_path_export
eval $ffeatools_pythonpath_export

#FFEA tools
echo "Install the FFEA python API? (y/n)"
read install_api
if [ "$install_api" == "y" ]; then
    ffeatools_api_pythonpath_export='export PYTHONPATH="$PYTHONPATH:'
    ffeatools_api_pythonpath_export+=$install_dir
    ffeatools_api_pythonpath_export+='/lib/python2.7/site-packages/"'
    eval $ffeatools_api_pythonpath_export
    echo $ffeatools_api_pythonpath_export >> ~/.bashrc
    cd ../ffea_src/ffeatools
    python setup.py install --prefix=$install_dir || exit 1
fi

#Phew! It's over.
echo "FFEA is now installed."
echo "To get started, visit our documentation at ffea.bitbucket.io/docs"
echo "To help us fund FFEA development, we humbly ask that you cite the research papers on the package."
