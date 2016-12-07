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
from math import ceil
import FFEA_traj

if len(sys.argv) != 6:
	sys.exit("Usage python " + os.path.basename(sys.argv[0]) + " [FFEA traj fname] [FFEA output traj fname] [frames to read] [First frame] [Last frame]")

traj_fname = sys.argv[1]
out_fname = sys.argv[2]
frames_to_read = int(sys.argv[3])
first_frame = int(sys.argv[4])
last_frame = int(sys.argv[5])

traj = FFEA_traj.FFEA_traj(traj_fname, frames_to_read, first_frame, last_frame, 1)
traj.write_traj_to_file(out_fname)
