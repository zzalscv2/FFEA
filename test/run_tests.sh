#!/bin/bash
# We're going to run some tests on all of the stuff in this folder, to make sure that ffea has installed correctly and is running and stuff!
scriptdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
ffeadir=$scriptdir/../bin

echo " "
echo "ffea: Usage..."
echo " "
sleep 3
$ffeadir/ffea

sleep 2

echo " "
echo "ffea: running with mass..."
echo " "
sleep 3
$ffeadir/ffea $scriptdir/sphere_63_120_mass/sphere_63_120_mass.ffea

echo " "
echo "ffea: running with no mass..."
echo " "
sleep 3
$ffeadir/ffea $scriptdir/sphere_63_120_nomass/sphere_63_120_nomass.ffea

echo " "
echo "ffea: running with multiple non-interacting blobs with no mass..."
echo " "
sleep 3
$ffeadir/ffea $scriptdir/sphere_63_120_two/sphere_63_120_two.ffea


echo " "
echo "ffea: running with multiple interacting blobs with mass..."
echo " "
sleep 3
$ffeadir/ffea $scriptdir/sphere_63_120_two_vdw/sphere_63_120_two_vdw.ffea
