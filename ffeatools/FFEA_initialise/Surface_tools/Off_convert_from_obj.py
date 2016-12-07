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
	sys.exit("Usage: python convert_obj_to_off.py [INPUT OBJ FNAME] [OUTPUT OFF FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line.split()[0] == "v":
		v.append(line.replace("v ",""))
	if line.split()[0] == "f":
		if "//" in line:
			line = line.replace("f ","").replace("//", " ")
			line = line.split()
			f.append(str(int(line[0])-1) + " " + str(int(line[2])-1) + " " + str(int(line[4])-1) + "\n")
		else:
			line = line.replace("f ", "")
			line = line.split()
			f.append(str(int(line[0])-1) + " " + str(int(line[1])-1) + " " + str(int(line[2])-1) + "\n")
obj.close()

off = open(sys.argv[2], "w")
off.write("OFF\n")
off.write(str(len(v)) + " " + str(len(f)) + " 0\n")
for vert in v:
	off.write(vert)
for face in f:
	off.write("3 " + face)
off.close()
