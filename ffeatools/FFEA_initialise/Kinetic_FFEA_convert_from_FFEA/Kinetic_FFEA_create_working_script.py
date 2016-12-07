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

if len(sys.argv) != 6:
	sys.exit("Usage: python " + sys.argv[0] + " [Original FFEA script fname] [Working script fname] [Working spring fname] [restart] [check]")

# Get args
py_script_dirname = os.path.dirname(sys.argv[0])
script_fname = sys.argv[1]
working_script_fname = sys.argv[2]
working_spring_fname = sys.argv[3]
restart = int(sys.argv[4])
check = int(sys.argv[5])

# Get lines from infile
fin = open(script_fname, "r")
lines = fin.readlines()
fin.close()

# Make new script
fout = open(working_script_fname, "w")

for line in lines:
	if "<spring>" in line or "</spring>" in line or "spring_fname" in line:
		continue
	elif "<interactions>" in line or "</interactions>" in line:
		continue
	elif "centroid_pos" in line:
		continue
	elif "</blob>" in line:
		#fout.write("\t\t<centroid_pos = (0.0,0.0,0.0)>\n")
		new_line = line
	elif "</system>" in line:
		break

	elif "=" not in line:
		new_line = line
	else:
		lvalue, rvalue = line.strip()[1:-1].split("=")
		if lvalue.strip() == "restart":
			if(restart == 0):
				new_line = line.replace(rvalue, " 0")	
			else:
				new_line = line.replace(rvalue, " 1")	
		elif lvalue.strip() == "check":
			new_line = line.replace(rvalue, " 100")
		elif lvalue.strip() == "num_steps":
			new_line = line.replace(rvalue, str(check * 100 * (restart + 1)))
		elif lvalue.strip() == "trajectory_out_fname" or lvalue.strip() == "measurement_out_fname":
			new_line = line.replace(rvalue, os.path.splitext(rvalue)[0] + "_working.out")
		elif lvalue.strip() == "calc_noise":
			new_line = line.replace(rvalue, " 0")
		elif lvalue.strip() == "calc_kinetics":
			new_line = line.replace(rvalue, " 0")
		elif lvalue.strip() == "stokes_visc":
			new_line = line.replace(rvalue, " 1e-04")
		elif lvalue.strip() == "pin":
			new_line = line.replace(rvalue, " working.pin")
			os.system("python " + py_script_dirname + "/../FFEA_pin_tools/pin_none.py working.pin")
		else:
			new_line = line
	fout.write(new_line)

if working_spring_fname != "none.spring":
	fout.write("\t<interactions>\n")
	fout.write("\t\t<springs>\n")
	fout.write("\t\t\t<spring_fname = " + working_spring_fname + ">\n")
	fout.write("\t\t</springs>\n")
	fout.write("\t</interactions>\n")

fout.write("</system>\n")
fout.close()
