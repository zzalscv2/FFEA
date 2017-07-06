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

import sys, os, subprocess
import __builtin__

import argparse as _argparse


# Set up argparse
parser = _argparse.ArgumentParser(description="Convert an FFEA trajectory to a pseudo-pdb system for PCA analysis")
parser.add_argument("i", help="Input PCZ file (.pcz)")
parser.add_argument("-n", action="store", nargs='?', default = '10', help="Number of Modes to Analyse")
parser.add_argument("-o", action="store", nargs='?', help="Output filename")

def FFEA_get_PCA_projections(infile, outfile, num_modes):

	# Check for problems
	base, ext = os.path.splitext(infile)

	if outfile == None:
		outfile = base + "_PCAproj"
	else:
		outfile = os.path.splitext(outfile)[0]

	outfiletemp = outfile + "_temp.dat"
	outfile = outfile + ".dat"

	if os.path.exists(outfile):
		print("Default output file ('" + outfile + "') already exists.\n")
		raise IOError

	try:
		num_modes = int(num_modes)
	except(ValueError):
		raise

	# Do some PCZ analysis
	# Print help to file and hack your way to num_evecs
	try:
		num_avail_modes = int(subprocess.check_output(["pyPczdump", "-i", infile, "-n"]).split("\n")[8][:-1].split()[-1])
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			raise OSError
		else:
			print("Unknown problem running 'pyPczdump. Perhaps conflicting versions (before and after 2.0)")
			raise IOError
	
	if num_modes > num_avail_modes:
		print("Too many modes requested. Defaulting to maximum (%d modes)" % (num_avail_modes))
		num_modes = num_avail_modes

	# Get some projections
	print("Writing projections to file...\n")
	fout = open(outfile, "w")
	for i in range(num_modes):
		proj = []
		sys.stdout.write("\r\tProjection " + str(i) + "...")
		try:
			subprocess.call(["pyPczdump", "-i", infile, "-p", str(i + 1), "-o", outfiletemp])
		except OSError as e:
			if e.errno == os.errno.ENOENT:
				raise OSError
			else:
				print("Unknown problem running 'pyPczdump. Perhaps conflicting versions (before and after 2.0)")
				raise IOError

		fin = open(outfiletemp, "r")
		lines = fin.readlines()
		fin.close()
		for l in lines:
			proj.append(float(l))

		for elem in proj:
			fout.write("%6.3f " % (elem))
		fout.write("\n")
	print("\tdone!")

	fout.close()
	print("done!")

	# Remove temp files
	subprocess.call(["rm", outfiletemp])

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    try:
	    args = parser.parse_args()
    except:
	somehelp = parser.format_help().split("\n", 1)[1]
	print somehelp
	sys.exit()

    try:
	FFEA_get_PCA_projections(args.i, args.o, args.n)
    except IOError:
        parser.print_help()
    except ValueError:
	print("'-n' must be an integer")
        parser.print_help()
    except TypeError:
        parser.print_help()
	print("\nLikely missing argument. Please try again :)\n")
    except OSError:
	print("\n'pyPczdump' program not found. Please add to your $PATH")
        parser.print_help()
