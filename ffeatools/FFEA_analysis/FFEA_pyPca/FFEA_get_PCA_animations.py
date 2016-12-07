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
import numpy as np
import FFEA_topology

if len(sys.argv) != 5:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .pcz file] [INPUT reference topology (_frame0.pdb)] [num_animations] [ffea scale]")

inpcz = sys.argv[1]
inref = sys.argv[2]
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
num_anim = int(sys.argv[3])
anim_basename = inpcz.split(".")[0]
ffea_scale = float(sys.argv[4])

# Get number of nodes for FFEA_stuff
print("Calculating Eigenvector Animations...")
for i in range(num_anim):
	anim_outfname = anim_basename + "_anim" + str(i) + ".pdb"
	anim_outfname_ffea = anim_basename + "_anim" + str(i) + ".ftj"
	print("\tEigenvector " + str(i) + ": Writing to " + os.path.splitext(os.path.basename(anim_outfname))[0] + ".pdb/.ftj ...")
	os.system("pyPczdump -i " + inpcz + " -m " + str(i) + " --pdb " + inref + " -o " + anim_outfname)
	os.system("python " + script_dir + "/../../FFEA_analysis/FFEA_traj_tools/PDB_convert_to_FFEA_trajectory.py " + anim_outfname + " " + anim_outfname_ffea + " " + str(ffea_scale))
	print("\tdone!")
print("done!")
