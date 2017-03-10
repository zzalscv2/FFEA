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
import FFEA_trajectory, FFEA_pdb, FFEA_kinetic_map

if len(sys.argv) != 6:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT FFEA trajectory (.ftj)] [INPUT Original .pdb] [INPUT FFEA map (.map)] [OUTPUT .pdb] [FFEA scale]")

# Get args
trajfname = sys.argv[1]
inpdbfname = sys.argv[2]
mapfname = sys.argv[3]
outpdbfname = sys.argv[4]
scale = float(sys.argv[5])

# Build objects
traj = FFEA_trajectory.FFEA_trajectory(trajfname)
pdb = FFEA_pdb.FFEA_pdb(inpdbfname)
kineticmap = FFEA_kinetic_map.FFEA_kinetic_map(mapfname)

print "Num_nodes in trajectory = ", traj.blob[0][0].num_nodes
print "Num_atoms in pdb = ", pdb.blob[0].num_atoms
print "Dimensions of map = ", kineticmap.num_rows, " x ", kineticmap.num_columns

# Make some empty frames
pdb.clear_position_data()
pdb.add_frames(traj.num_frames)


# Apply map to each frame
print("Applying Map...")
for i in range(traj.num_frames):
	sys.stdout.write("\r%d%% read" % (int(i * 100.0 / traj.num_frames)))
	sys.stdout.flush()
	new_pos = (1.0 / scale) * kineticmap.apply_sparse(traj.blob[0][0].frame[i])
	pdb.blob[0].set_pos(i, new_pos)

print("\ndone!")
print("Wrinting PDB to " + outpdbfname + "...")
pdb.write_to_file(outpdbfname)
print("done!")
