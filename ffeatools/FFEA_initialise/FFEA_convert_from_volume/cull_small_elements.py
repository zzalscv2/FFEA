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

def jac(n1, n2, n3, n4, scale):
	return [[n2[0] * scale - n1[0] * scale, n2[1] * scale - n1[1] * scale, n2[2] * scale - n1[2] * scale], [n3[0] * scale - n1[0] * scale, n3[1] * scale - n1[1] * scale, n3[2] * scale - n1[2] * scale], [n4[0] * scale - n1[0] * scale, n4[1] * scale - n1[1] * scale, n4[2] * scale - n1[2] * scale]]

def get_vol_from_J(J):
	det = J[0][0] * (J[2][2]*J[1][1] - J[2][1]*J[1][2]) + J[1][0] * (J[2][1]*J[0][2] - J[2][2]*J[0][1]) + J[2][0] * (J[1][2]*J[0][1] - J[1][1]*J[0][2]);
	return abs(det)/6.0

def vol(element, nodes, scale):
	J = jac(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]], scale)
	return get_vol_from_J(J)

def is_i_in_use(i, elements):
	for el in elements:
		for n in el:
			if i == n:
				return True
	return False

if len(sys.argv) != 5:
	sys.exit("Usage: python cull_small_elements.py [.NODE FILE] [.ELE FILE] [CULL THRESHOLD VOL] [SCALE]")

inputnode = open(sys.argv[1],"r")
inputelem = open(sys.argv[2],"r")
cull_thresh = float(sys.argv[3])
scale = float(sys.argv[4])

# Get number of nodes
line = inputnode.readline().split()
num_nodes = int(line[0])
print "num_nodes = ", num_nodes

nodes = []
# Get nodes
for i in range(num_nodes):
	line = inputnode.readline().split()
	nodes.append([float(line[1]), float(line[2]), float(line[3])])

inputnode.close()

# Get number of elements
line = inputelem.readline().split()
num_elem = int(line[0])
print "num_elem = ", num_elem

# Get elements, culling those below the specified volume
elements = []
for i in range(num_elem):
	line = inputelem.readline().split()
	element = [int(line[1]), int(line[2]), int(line[3]), int(line[4])]
	if vol(element, nodes, scale) > cull_thresh:
		elements.append(element)

inputelem.close()

new_num_elem = len(elements)
print "Culled " + str(num_elem - new_num_elem) + " elements."

# Get how many nodes are now unused
for i in range(num_nodes - 1, -1, -1):
	if is_i_in_use(i, elements) == False:
		print "node " + str(i) + " is now unused"
		# remove node from node list
		nodes.pop(i)
		# renumber indices in elements
		for j in range(new_num_elem):
			for k in range(4):
				if (elements[j])[k] > i:
					(elements[j])[k] -= 1

new_num_nodes = len(nodes)
print "Culled " + str(num_nodes - new_num_nodes) + " nodes."

outputnode = open(sys.argv[1], "w")
outputnode.write(str(new_num_nodes) + " 3 0 0\n")
i = 0
for n in nodes:
	outputnode.write(str(i) + " " + str(n[0]) + " " + str(n[1]) + " " + str(n[2]) + "\n")
	i += 1
outputnode.close()

outputelem = open(sys.argv[2], "w")
outputelem.write(str(new_num_elem) + " 4 0\n")
i = 0
for el in elements:
	outputelem.write(str(i) + " " + str(el[0]) + " " + str(el[1]) + " " + str(el[2]) + " " + str(el[3]) + "\n")
	i += 1
outputelem.close()

print "cull_small_elements.py -> Done."
