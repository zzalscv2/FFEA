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
	sys.exit("Usage: python " + sys.argv[0] + "[INPUT input_fname (.obj)] [OUTPUT output_fname (.gts)]\n")
	
ifile = sys.argv[1]
ofile = sys.argv[2]
verts = 0
edges = 0
faces = 0

sta = open(ifile,'r')
stb = open(ofile,'w')
for i in sta: 
  if i[:2] == "v ":
    verts += 1 
  elif i[:2] == "f ":
    faces += 1  
    edges += 3

# go back to the top,
sta.seek(0)
#   and write the vertices:
stb.write(str(verts) + ' ' + str(edges) + ' ' + str(faces) + ' GtsSurface GtsFace GtsEdge GtsVertex\n')
for i in sta: 
  if i[:2] == "v ":
    stb.write(i[2:])

# go back to the top,
sta.seek(0)
#   and write triangular stuff:
for i in sta:
  if i[:2] == "f ":
    l = i.split()
    v1 = str(int(l[1].split("//")[0]))
    v2 = str(int(l[2].split("//")[0]))
    v3 = str(int(l[3].split("//")[0]))
    stb.write(v1 + " " + v2 + "\n")
    stb.write(v2 + " " + v3 + "\n")
    stb.write(v3 + " " + v1 + "\n")

sta.close()

##### and finally: 
for i in range(faces):
  e1 = 3*i + 1 
  stb.write(str(e1) + " " + str(e1 + 1) + " " + str(e1 + 2) + "\n")

stb.close()
