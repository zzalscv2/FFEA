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

import sys, os
import FFEA_script

if len(sys.argv) != 6:
	sys.exit("Usage: python FFEA_generate_kinetic_maps.py [INPUT .ffea file] [Blob index 1] [Conformation Index 1] [Blob index 2] [Conformation Index 2]\n")

# Get args
infname = sys.argv[1]
bindex = [int(sys.argv[2]), int(sys.argv[4])]
cindex = [int(sys.argv[3]), int(sys.argv[5])]

scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
ffeadir = os.path.dirname(os.path.abspath(infname))

# Get a script and the needed files from it
script = FFEA_script.FFEA_script(infname)

# Make the structures overlap
os.system("python " + scriptdir + "/FFEA_make_structures_overlap.py " + infname + " " + str(bindex[0]) + " " + str(cindex[0]) + " " + str(bindex[1]) + " " + str(cindex[1]))

# Make the maps!
basenode = ffeadir + "/blob0_overlap.node"
targetnode = ffeadir + "/blob1_overlap.node"
basetop = script.blob[bindex[0]].conformation[cindex[0]].topology
targettop = script.blob[bindex[1]].conformation[cindex[1]].topology

basetotargetmap = infname.split(".")[0] + "_b" + str(bindex[0]) + "c" + str(cindex[0]) + "to" + "b" + str(bindex[1]) + "c" + str(cindex[1]) + ".map"
targettobasemap = infname.split(".")[0] + "_b" + str(bindex[1]) + "c" + str(cindex[1]) + "to" + "b" + str(bindex[0]) + "c" + str(cindex[0]) + ".map"

#basetotargetsparsemap = basetotargetmap.split(".")[0] + "_sparse.map"
#targettobasesparsemap = targettobasemap.split(".")[0] + "_sparse.map"

os.system(scriptdir + "/make_structure_map -i " + basenode + " -t " + basetop + " -o " + targetnode + " -m " + basetotargetmap)
os.system(scriptdir + "/make_structure_map -i " + targetnode + " -t " + targettop + " -o " + basenode + " -m " + targettobasemap)

# Make them sparse
#os.system("python " + scriptdir + "/FFEA_convert_kinetic_map_to_sparse.py " + basetotargetmap + " " + basetotargetsparsemap)
#os.system("python " + scriptdir + "/FFEA_convert_kinetic_map_to_sparse.py " + targettobasemap + " " + targettobasesparsemap)

# Remove the unneeded stuff
#os.system("rm " + basenode + " " + targetnode + " " + basetotargetmap + " " + targettobasemap)
#os.system("rm " + basetotargetmap + " " + targettobasemap)
#os.system("rm " + basenode + " " + targetnode)

# Finalise
print("Maps created:")
print("\t" + basetotargetmap)
print("\t" + targettobasemap)
