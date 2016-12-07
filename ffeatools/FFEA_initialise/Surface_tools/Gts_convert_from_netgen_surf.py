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
if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + "[INPUT input_fname (.surf)] [OUTPUT output_fname (.gts)]\n")
	
ifile = sys.argv[1]
ofile = sys.argv[2]
verts = 0
edges = 0
faces = 0

sta = open(ifile,'r')
stb = open(ofile,'w')
sta.readline()
verts = int(sta.readline())
for i in range(verts):
	sta.readline()
faces = int(sta.readline())
edges = faces * 3

# go back to the top,
sta.seek(0)
#   and write the vertices:
stb.write(str(verts) + ' ' + str(edges) + ' ' + str(faces) + ' GtsSurface GtsFace GtsEdge GtsVertex\n')
sta.readline()
sta.readline()
for i in range(verts): 
	stb.write(sta.readline())

#   and write triangular stuff:
sta.readline()
for i in range(faces):
	l = sta.readline().split()
	v1 = l[0].strip()
	v2 = l[1].strip()
	v3 = l[2].strip()
	stb.write(v1 + " " + v2 + "\n")
	stb.write(v2 + " " + v3 + "\n")
	stb.write(v3 + " " + v1 + "\n")

sta.close()

##### and finally: 
for i in range(faces):
  e1 = 3*i + 1 
  stb.write(str(e1) + " " + str(e1 + 1) + " " + str(e1 + 2) + "\n")

stb.close()
