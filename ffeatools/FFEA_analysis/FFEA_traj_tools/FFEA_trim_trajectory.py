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
import FFEA_trajectory
import argparse as _argparse
import __builtin__

parser = _argparse.ArgumentParser(description="FFEA Trim Trajectory")
parser.add_argument("traj_fname", action="store", help="Input trajectory file (.ftj)")
parser.add_argument("out_fname", action="store", help="Output trajectory file (.ftj)")
parser.add_argument("frames_to_read", action="store", type=int, help="Number of frames to read")
parser.add_argument("trim_percent", action="store", type=float, help="Percentage to keep")

def trim_trajectory(traj_fname, frames_to_read, trim_percent):
    """
    Remove frames from start of an FFEA_trajectory file.
    In:
    traj_fname an ffea trajectory file (note: a file, not an instance!)
    out_fname: an output filename
    frames_to_read: the number of frames to read,
    thin_percent: percentage of the file to keep
    Retruns:
    an FFEA trajectory object.
    """
    
    if trim_percent < 0 or trim_percent > 100:
        sys.exit("Error. Percentage must be between 0 and 100. You used %f\n" % (trim_percent))
    
    if trim_percent < 1:
        verify = raw_input("Percentage to keep was %f. Did you mean %f (y/n)?" % (trim_percent, thin_percent * 100))
        if verify.lower() == "y":
            trim_percent *= 100
    
    trim_percent /= 100.0

    total_num_frames = FFEA_trajectory.get_num_frames(traj_fname)
    if frames_to_read > total_num_frames:
	frames_to_read = total_num_frames

    traj = FFEA_trajectory.FFEA_trajectory(traj_fname, start = frames_to_read * (1 - trim_percent), num_frames_to_read = frames_to_read)
    return traj
    
if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    # Get args and build objects
    print args.traj_fname
    if not os.path.exists(args.traj_fname):
        raise IOError("Trajectory file specified doesnae exist.")
    out_traj = trim_trajectory(args.traj_fname, args.frames_to_read, args.trim_percent, )
    out_traj.write_to_file(args.out_fname)
