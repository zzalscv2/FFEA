#!/bin/bash
# We're going to run some tests on all of the stuff in this folder, to make sure that ffea has installed correctly and is running and stuff!
scriptdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
ffeadir=$scriptdir/../../src
export OMP_NUM_THREADS=4

echo " "
echo "ffea: Usage..."
echo " "
sleep 3
ffea

echo " "
echo "ffea: Simulation of Flat Cuboid..."
echo "          num_elements = 179"
echo "          num_nodes    = 452"
echo " "
sleep 3
cd $scriptdir/cuboid_flat
ffea cuboid_flat.ffea
cd - 

echo " "
echo "ffea: Simulation of Multiple Flat Cuboids Interacting..."
echo "          num_elements = 179"
echo "          num_nodes    = 1452"
echo " "
sleep 3
cd $scriptdir/cuboid_flat/
ffea cuboid_flat_multiple.ffea
cd - 

echo " "
echo "ffea: running with mass..."
echo " "
sleep 3
cd $scriptdir/sphere_63_120_mass/
ffea sphere_63_120_mass.ffea
cd -

echo " "
echo "ffea: running with no mass..."
echo " "
sleep 3
cd $scriptdir/sphere_63_120_nomass/
ffea sphere_63_120_nomass.ffea
cd -

echo " "
echo "ffea: running with multiple non-interacting blobs with no mass..."
echo " "
sleep 3
cd $scriptdir/sphere_63_120_two/
ffea sphere_63_120_two.ffea
cd - 

echo " "
echo "ffea: running with multiple interacting blobs with mass..."
echo " "
sleep 3
cd $scriptdir/sphere_63_120_two_vdw/
ffea sphere_63_120_two_vdw.ffea
cd -

echo " "
echo "ffea: running elastic network model of flat cuboid..."
echo " "
sleep 3
cd $scriptdir/cuboid_flat/
ffea -m 1 cuboid_flat_enm.ffea << EOF
0
q
5
EOF
cd -

echo " "
echo "ffea: running dynamic mode model of flat cuboid..."
echo " "
sleep 3
cd $scriptdir/cuboid_flat/
ffea -m 2 cuboid_flat_dmm.ffea << EOF
0
q
5
EOF
cd -
