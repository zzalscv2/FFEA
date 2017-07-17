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
import FFEA_script, FFEA_trajectory, FFEA_pdb
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Convert an FFEA trajectory to a pseudo-pdb system for PCA analysis")
parser.add_argument("i", help="Input file (.ffea / .ftj)")
parser.add_argument("-o", action="store", nargs='?', help="Output file (.pdb).")
parser.add_argument("-ind", action="store", nargs='?', default = '0', help="Blob index")
parser.add_argument("-frames", action="store", nargs='?', default = None, help="Number of frames to read")

def FFEA_convert_to_PCAsystem(infile, outfile, frames, bindex):

	# Check for problems
	base, ext = os.path.splitext(infile)

	if outfile == None:
		outfile = base + "_PCAready.pdb"
		if os.path.exists(outfile):
			print("Default output file '" + outfile + "' already exists.\n")
			raise IOError

	# Get an input file
	try:
		if ext == ".ffea":
			script = FFEA_script.FFEA_script(infile)
			if frames == None:
				traj = script.load_trajectory()
			else:
				traj = script.load_trajectory(int(frames))

		elif ext == ".ftj":
			if frames == None:
				traj = FFEA_trajectory.FFEA_trajectory(infile)
			else:
				traj = FFEA_trajectory.FFEA_trajectory(infile, num_frames_to_read = int(frames))
	except:
		print("Could not read from input file.")
		raise IOError
	
	# Make traj into a single blob traj
	while(True):
		try:
			if traj.blob[int(bindex)][0].motion_state != "DYNAMIC":
				print("Selected a STATIC blob. Unable to process.")
				raise IndexError

			traj.set_single_blob(int(bindex))
			break

		except(IndexError):
			print("\nIncorrect blob index chosen.")
			try:
				bindex = input("Please enter a new blob index (less than " + str(traj.num_blobs) + "):")
			except:
				print("Could not read blob index")
				raise IndexError

	# Make the thing into a pdb
	pdb = FFEA_pdb.FFEA_pdb("")
	pdb.build_from_traj(traj)

	# Write first frame to a file (for a PCA topology file)
	base, ext = os.path.splitext(outfile)
	outfilef0 = base + "_frame0" + ext
	pdb.write_to_file(outfilef0, frames = [0,1])

	# Now write whole lot to a file
	pdb.write_to_file(outfile)

	base, ext = os.path.splitext(infile)
	if ext == ".ffea":
		print("\nScript Filename - " + infile)
		print("Original Traj Filename - " + script.params.trajectory_out_fname)
	else:
		print("Original Traj Filename - " + infile)	

	print("NEW first frame pdb (send to pcazip as topology)- " + outfilef0)
	print("NEW pdb trajectory (send to pcazip as trajectory)- " + outfile)
	print("\n")

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    try:
	    args = parser.parse_args()
    except:
	somehelp = parser.format_help().split("\n", 1)[1]
	print somehelp
	sys.exit()

    try:
        FFEA_convert_to_PCAsystem(args.i, args.o, args.frames, args.ind)
    except IOError:
        parser.print_help()
    except TypeError:
        parser.print_help()
	print("\nLikely missing argument. Please try again :)\n")
