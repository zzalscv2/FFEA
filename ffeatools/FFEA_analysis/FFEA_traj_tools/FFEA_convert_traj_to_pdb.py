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
import FFEA_trajectory, FFEA_pdb

if len(sys.argv) != 5:
	sys.exit("Usage: python FFEA_convert_traj_to_pdb.py [INPUT FFEA traj fname (.out)] [OUTPUT .pdb fname] [num_frames_to_read] [FFEA scale]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
num_frames_to_read = int(sys.argv[3])
scale = 1.0 / float(sys.argv[4])	# Invert as pdb is in angstroms

# Build objects
traj = FFEA_trajectory.FFEA_trajectory(infname)
pdb = FFEA_pdb.FFEA_pdb("")
pdb.build_from_traj(traj, scale = scale)
pdb.write_to_file(outfname)
