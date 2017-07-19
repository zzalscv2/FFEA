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

# Set up argparse
parser = _argparse.ArgumentParser(description="Extract a subsection of an FFEA trajectory")
parser.add_argument("i", help="Input trajectory file (.ftj)")
parser.add_argument("-s", action="store", nargs=1, default = '0', help="Snapshot to begin from")
parser.add_argument("-f", action="store", nargs=1, help="Snapshot to end from")
parser.add_argument("-o", action="store", nargs='?', help="Output filename")

def FFEA_split_trajectory(infile, outfile, start, end):

	# Check for problems
	base, ext = os.path.splitext(infile)

	if end == None:
		print("Error: We need a final frame index (-f) to define a subsection of the trajectory.")
		raise ValueError

	start = int(start[0])
	end = int(end[0])

	if start >= end:
		print("Error: start value greater than end value. Please reenter.")
		raise ValueError

	# Read only the amount of trajectory that we were asked for
	try:
		traj = FFEA_trajectory.FFEA_trajectory(infile, num_frames_to_read = end - start, start = start)
	except:
		raise

	if end - start != traj.num_frames:
		if sys.version_info[0] < 3:
			inputter = raw_input
		else:
			inputter = input

		ans = inputter("Number of frames within trajectory, %d, less than specified range, %d - %d. Is this ok (y/n)?: " % (traj.num_frames, start, end))
	
		try:
			if str(ans).lower() == "n":
				print("Will not return a new trajectory file. Please choose a new range next time :)")
			else:
				raise ValueError
		except:
			print("Assuming 'y' was entered. Get ready for a new trajectory file!")
			end = start + traj.num_frames

	if outfile == None:
		outfile = base + "_extracted" + str(start) + "-" + str(end) + ".ftj"

	traj.write_to_file(outfile)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
	try:
		args = parser.parse_args()
	except:
		somehelp = parser.format_help().split("\n", 1)[1]
		print(somehelp)
		sys.exit()

	try:
		FFEA_split_trajectory(args.i, args.o, args.s, args.f)
	except IOError:
	        parser.print_help()
	except ValueError:
		parser.print_help()
