#!/bin/bash
# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

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
