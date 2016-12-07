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

import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT OBJ FNAME] [OUTPUT SURF FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line[0] != "v" and line[0] != "f":
		continue
	one_line = line.split()
	if one_line[0] == "v":
		v.append(line.replace("v ",""))
	if one_line[0] == "f":
		index = [0,0,0]
		index[0] = int(one_line[1].split("//")[0])
		index[1] = int(one_line[2].split("//")[0])
		index[2] = int(one_line[3].split("//")[0])
		for i in range(len(index)):
			if index[i] < 0:
				index[i] = -1 * index[i]

		face_str = str(index[0]) + " " + str(index[1]) + " " + str(index[2]) + "\n"
		f.append(face_str)
obj.close()

surf = open(sys.argv[2], "w")
surf.write("surfacemesh\n")
surf.write(str(len(v)) + "\n")
for vert in v:
	surf.write(vert)
surf.write(str(len(f)) + "\n")
for face in f:
	surf.write(face)
surf.close()
