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
import FFEA_script, FFEA_trajectory, FFEA_measurement
import argparse as _argparse
import __builtin__

parser = _argparse.ArgumentParser(description="FFEA Thin Trajectory")
parser.add_argument("script_fname", action="store", help="Input script file (.ffea)")
parser.add_argument("out_fname", action="store", help="Output trajectory/measurement file (.ftj/.fm/.fdm)")
parser.add_argument("frames_to_read", action="store", type=int, help="Number of frames to read")
parser.add_argument("thin_percent", action="store", type=float, help="Percentage to keep")

def thin_system(script_fname, frames_to_read, thin_percent):
    """
    Remove frames from an FFEA_trajectory file.
    In:
    traj_fname an ffea trajectory file (note: a file, not an instance!)
    out_fname: an output filename
    frames_to_read: the number of frames to read,
    thin_percent: percentage of the file to keep
    Retruns:
    an FFEA trajectory object.
    """
    
    if thin_percent < 0 or thin_percent > 100:
        sys.exit("Error. Percentage must be between 0 and 100. You used %f\n" % (thin_percent))
    
    if thin_percent < 1:
        verify = raw_input("Percentage to keep was %f. Did you mean %f (y/n)?" % (thin_percent, thin_percent * 100))
        if verify.lower() == "y":
            thin_percent *= 100
    
    frame_rate = ceil(100 / thin_percent)
    script = FFEA_script.FFEA_script(script_fname)
    traj = FFEA_trajectory.FFEA_trajectory(script.params.trajectory_out_fname, frame_rate = frame_rate, num_frames_to_read = frames_to_read)
    meas = FFEA_measurement.FFEA_measurement(script.params.measurement_out_fname, frame_rate = frame_rate, num_frames_to_read = frames_to_read)
    return traj, meas
    
if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    # Get args and build objects
    if not os.path.exists(args.script_fname):
        raise IOError("Script file specified doesn't exist.")
    out_traj, out_meas = thin_system(args.script_fname, args.frames_to_read, args.thin_percent, )

    # Get output_fnames and write out
    out_bfname = os.path.splitext(args.out_fname)[0]
    out_tfname = out_bfname + ".ftj"
    out_mfname = out_bfname + ".fm"
    out_traj.write_to_file(out_tfname)
    out_meas.write_to_file(out_mfname)
