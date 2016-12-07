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

import sys

if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT OFF FILE] [OUTPUT NETGEN SURFACE FILE]")

file_in = open(sys.argv[1],"r")
file_out = open(sys.argv[2],"w")

# read in the OFF line
line = file_in.readline()
if "OFF" not in line:
	print "Not an OFF file (first line is not OFF)"

# write the header line in the output file
file_out.write("surfacemesh\n")

# read in the line with the number of nodes, number of faces, and some mysterious number I have no idea about
line = file_in.readline()
a = line.split()
num_nodes = int(a[0])
num_faces = int(a[1])
print "num_nodes = " + str(num_nodes)
print "num_faces = " + str(num_faces)

# write the nodes section of the output file
file_out.write(str(num_nodes) + "\n")
for i in range(num_nodes):
	line = file_in.readline()
	file_out.write(line)

# write the faces section of the output file
file_out.write(str(num_faces) + "\n")
for i in range(num_faces):
	line = file_in.readline()
	a = line.split()
	n1 = str(int(a[1]) + 1)
	n2 = str(int(a[2]) + 1)
	n3 = str(int(a[3]) + 1)

	# miss out the mysterious first number
	file_out.write(n1 + " " + n2 + " " + n3 + "\n")


file_in.close()
file_out.close()
