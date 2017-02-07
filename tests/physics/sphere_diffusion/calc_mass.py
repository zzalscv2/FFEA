import sys, os
import FFEA_script
import numpy as np
try:
	script = FFEA_script.FFEA_script(sys.argv[1])
except:
	sys.exit("Enter a script fbname please")

try:
	top = script.load_topology(0)
	mat = script.load_material(0)
	node = script.load_node(0)
except:
	sys.exit("Indexing starts from zero!")

print top.calc_mass(mat=mat, node=node, scale = script.blob[0].scale)

matlin = np.array([[0.0 for i in range(4)]for j in range(4)])
matquad = np.array([[0.0 for i in range(12)]for j in range(12)])
mass = 1

for i in range(12):
	for j in range(12):
		if i < 4 and j < 4:
			if i == j:
				matlin[i][j] = 0.1 * mass
				matquad[i][j] = 0.1 * mass
			else:
				matlin[i][j] = 0.05 * mass
				matquad[i][j] = 0.05 * mass
		else:
			if i == j:
				matquad[i][j] = 1.0

print matquad
evals, evecs = np.linalg.eig(matquad)

print evals
print evecs

print np.dot(matlin, np.array([1,1,1,1]))

#a = np.array([[1,0,0],[0,2,0],[0,0,3]])
#evals, evecs = np.linalg.eig(a)

#print evals
#print evecs

#print np.dot(a, evecs[2])
