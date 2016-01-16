#!/bin/bash
# We're going to run some tests on all of the stuff in this folder, to make sure that ffea has installed correctly and is running and stuff!
scriptdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
ffeadir=$scriptdir/../build/src
export OMP_NUM_THREADS=4

echo " "
echo "ffea: Usage..."
echo " "
sleep 3
$ffeadir/ffea

echo " "
echo "ffea: Simulation of Flat Cuboid..."
echo "          num_elements = 465"
echo "          num_nodes    = 1025"
echo " "
sleep 3
$ffeadir/ffea/cuboid_flat/scripts/cuboid_flat.ffea

echo " "
echo "ffea: Simulation of Multiple Flat Cuboids Interacting..."
echo "          num_elements = 465"
echo "          num_nodes    = 1025"
echo " "
sleep 3
$ffeadir/ffea/cuboid_flat/scripts/cuboid_flat_multiple.ffea

#echo " "
#echo "ffea: running with mass..."
#echo " "
#sleep 3
#$ffeadir/ffea $scriptdir/sphere_63_120_mass/sphere_63_120_mass.ffea

#echo " "
#echo "ffea: running with no mass..."
#echo " "
#sleep 3
#$ffeadir/ffea $scriptdir/sphere_63_120_nomass/sphere_63_120_nomass.ffea

#echo " "
#echo "ffea: running with multiple non-interacting blobs with no mass..."
#echo " "
#sleep 3
#$ffeadir/ffea $scriptdir/sphere_63_120_two/sphere_63_120_two.ffea


#echo " "
#echo "ffea: running with multiple interacting blobs with mass..."
#echo " "
#sleep 3
#$ffeadir/ffea $scriptdir/sphere_63_120_two_vdw/sphere_63_120_two_vdw.ffea

#echo " "
#echo "ffea: running elastic network model of flat cuboid..."
#echo " "
#sleep 3
#$ffeadir/ffea -m 1 $scriptdir/cuboid_flat/scripts/cuboid_flat_enm.ffea << EOF
#0
#q
#5
#EOF

#echo " "
#echo "ffea: running elastic network model of flat cuboid..."
#echo " "
#sleep 3
#$ffeadir/ffea -m 2 $scriptdir/cuboid_flat/scripts/cuboid_flat_dmm.ffea << EOF
#0
#q
#5
#EOF
