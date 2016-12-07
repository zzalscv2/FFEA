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

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT OBJ FNAME] [OUTPUT FFEA SURF FNAME] [OUTPUT FFEA NODES FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line == "\n":
		continue
	if line.split()[0] == "v":
		v.append(line.replace("v ",""))
	if line.split()[0] == "f":
		f.append(line.replace("f ",""))
obj.close()

surf = open(sys.argv[2], "w")
surf.write("FFEA surface file\n")
surf.write("num_surface_faces " + str(len(f)) + "\n")
surf.write("faces:\n")
for face in f:
	face_split = face.split()
	index = [0,0,0]
	index[0] = face_split[0].split("//")[0]
	index[1] = face_split[1].split("//")[0]
	index[2] = face_split[2].split("//")[0]
	face_str = "0 " + str(int(index[0]) - 1) + " " + str(int(index[1]) - 1) + " " + str(int(index[2]) - 1) + "\n"
	surf.write(face_str)
surf.close()

node = open(sys.argv[3], "w")
node.write("FFEA node file\n")
num_nodes = len(v)
node.write("num_nodes " + str(num_nodes) + "\n")
node.write("num_surface_nodes " + str(num_nodes) + "\n")
node.write("num_interior_nodes 0\n")
node.write("surface nodes:\n")
for vert in v:
	node.write(vert)
node.write("interior nodes:\n")
node.close()
