import sys, math

class Node:
	def __init__(self):
		self.x = 0.0
		self.y = 0.0
		self.z = 0.0

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def mag(self):
		return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

# Importing functions
def import_nodes_from_ffea(nodes_fname):

	ffea_nodes = []
	infile = open(nodes_fname, "r")
	line = infile.readline().strip()
	if line != "walrus node file" and line != "ffea node file":
		sys.exit(rvalue + " is not an ffea node file, or is formatted incorrectly\n") 

	num_nodes = int(infile.readline().split()[1])
	num_surface_nodes = int(infile.readline().split()[1])
	num_interior_nodes = int(infile.readline().split()[1])
	infile.readline()

	for i in range(num_surface_nodes):
		sline = infile.readline().split()
		ffea_nodes.append(Node(float(sline[0]), float(sline[1]), float(sline[2])))
	
	infile.readline()
	for i in range(num_interior_nodes):
		sline = infile.readline().split()
		ffea_nodes.append(Node(float(sline[0]), float(sline[1]), float(sline[2])))

	infile.close()
	return ffea_nodes
	
		
