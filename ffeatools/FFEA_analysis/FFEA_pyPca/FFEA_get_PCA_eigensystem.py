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
import FFEA_topology
import numpy as np
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Convert an FFEA trajectory to a pseudo-pdb system for PCA analysis")
parser.add_argument("i", help="Input PCZ file (.pcz)")
parser.add_argument("t", help="Input FFEA topology file (.top)")
parser.add_argument("-n", action="store", nargs='?', default = '10', help="Number of Modes to Analyse")
parser.add_argument("-o", action="store", nargs='?', help="Output filename")

def FFEA_get_PCA_eigensystem(infile, topfile, outfile, num_modes):

	# Check for problems
	base, ext = os.path.splitext(infile)

	if outfile == None:
		outfile = base + "_PCAeigensys"
	else:
		outfile = os.path.splitext(outfile)[0]

	outfilevec = outfile + ".evecs"
	outfileval = outfile + ".evals"
	tempoutfilevec = outfile + "_temp.evecs"
	if os.path.exists(outfilevec) or os.path.exists(outfilevec):
		print("Default output file ('" + outfilevec + "') or ('" + outfileval + "') already exists.\n")
		raise IOError

	try:
		num_modes = int(num_modes)
	except(ValueError):
		raise

	# Read topology file
	try:
		top = FFEA_topology.FFEA_topology(topfile)
	except:
		print("Could not read topology file.")
		raise IOError
	lin = top.get_linear_nodes()

	# Do some PCZ analysis

	# Check version (for some reason, it's written to stderr :/)
	p = subprocess.Popen(["pyPczdump", "--version"], stderr=subprocess.PIPE)
	sys.stderr.flush()
	pyPczver = p.communicate()[1].strip()
	sys.stdout.write("Found pyPczdump version " + pyPczver + "\n\n")
	pyPczver = [int(bit) for bit in pyPczver.split(".")]

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

	# Evals
	print("Getting Eigenvalues...")
	try:
		subprocess.call(["pyPczdump", "-i", infile, "-l", "-o", outfileval])
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			raise OSError
		else:
			print("Unknown problem running 'pyPczdump. Perhaps conflicting versions (before and after 2.0)")
			raise IOError

	# Reduce to number of evals requested
	fin = open(outfileval, "r")
	lines = fin.readlines()
	fin.close()
	
	fin = open(outfileval, "w")
	for i in range(num_modes):
		fin.write(lines[i])
	fin.close()
	print("...done!\n")

	# Evecs
	print("Getting Eigenvectors...\n")
	fout = open(outfilevec, "w")
	for i in range(int(num_modes)):
		sys.stdout.write("\rEigenvector %d" % (i + 1))
		sys.stdout.flush()
		evec = []
		if(pyPczver[0] >= 2):
			try:
				subprocess.call(["pyPczdump", "-i", infile, "-e", str(i + 1), "-o", tempoutfilevec])
			except OSError as e:
				if e.errno == os.errno.ENOENT:
					raise OSError
				else:
					print("Unknown problem running 'pyPczdump. Perhaps conflicting versions (before and after 2.0)")
					raise IOError
		else:
			try:
				subprocess.call(["pyPczdump", "-i", infile, "-e", str(i), "-o", tempoutfilevec])
			except OSError as e:
				if e.errno == os.errno.ENOENT:
					raise OSError
				else:
					print("Unknown problem running 'pyPczdump. Perhaps conflicting versions (before and after 2.0)")
					raise IOError

		with open(tempoutfilevec, "r") as fin:
			lines = fin.readlines()

		num_nodes = len(lines) / 3
		for j in range(num_nodes):
			for k in range(3):
				if j in lin:
					evec.append(float(lines[3 * j + k]))			

		evec = np.array(evec)
		evec /= np.linalg.norm(evec)

		for elem in evec:
			fout.write("%6.3f " % (elem))
		fout.write("\n")

	print("\n...done!")
	fout.close()

	# Remove temp files
	subprocess.call(["rm", tempoutfilevec])

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    try:
        args = parser.parse_args()
    except:
        somehelp = parser.format_help().split("\n", 1)[1]
        print(somehelp)
        sys.exit()

    try:
        FFEA_get_PCA_eigensystem(args.i, args.t, args.o, args.n)
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
