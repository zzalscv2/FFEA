import sys, os
from math import ceil
import FFEA_trajectory
import argparse as _argparse

parser = _argparse.ArgumentParser(description="FFEA Thin Trajectory")
parser.add_argument("traj_fname", action="store", help="Input trajectory file (.ftj)")
parser.add_argument("out_fname", action="store", help="Output trajectory file (.ftj)")
parser.add_argument("frames_to_read", action="store", type=int, help="Number of frames to read")
parser.add_argument("thin_percent", action="store", type=float, help="Percentage to keep")

def thin_trajectory(traj_fname, frames_to_read, thin_percent):
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
    traj = FFEA_trajectory.FFEA_trajectory(traj_fname, frame_rate = frame_rate, num_frames_to_read = frames_to_read)
    return traj
    
if sys.stdin.isatty():
    args = parser.parse_args()
    # Get args and build objects
    print args.traj_fname
    if not os.path.exists(args.traj_fname):
        raise IOError("Trajectory file specified doesnae exist.")
    out_traj = thin_trajectory(args.traj_fname, args.frames_to_read, args.thin_percent, )
    out_traj.write_to_file(args.out_fname)
