import sys
from Vectors import vector3

class blob():
	def __init__(self):
		self.state = " "
		self.nodes = " "
		self.vol = " "
		self.scale = 0.0

class params():
	def __init__(self):
		self.es_h = 0
		self.kappa = 0.0
		self.es_N_x = 0
		self.es_N_y = 0
		self.es_N_z = 0
		self.num_blobs = 0
class vector3D():
	def __init__(self):
		self.x = 0.0
		self.y = 0.0
		self.z = 0.0

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		
def read_param_block(lines, myparams):
	for line in lines:
		if line == "<param>\n":
			param_index = lines.index("<param>\n")
			break
	
	check = 0
	while line != "</param>\n":
		param_index = param_index + 1
		if line != "\n":
			line = lines[param_index]
			sline = line.split()
			if sline[0][1:] == "es_h":
				myparams.es_h = sline[2].split(">")[0]
				check = check + 1
			elif sline[0][1:] == "kappa":
				myparams.kappa = sline[2].split(">")[0]	
				check = check + 1
			elif sline[0][1:] == "es_N_x":
				myparams.es_N_x = sline[2].split(">")[0]
				check = check + 1
			elif sline[0][1:] == "es_N_y":
				myparams.es_N_y = sline[2].split(">")[0]
				check = check + 1
			elif sline[0][1:] == "es_N_z":
				myparams.es_N_z = sline[2].split(">")[0]
				check = check + 1
			elif sline[0][1:] == "num_blobs":
				myparams.num_blobs = sline[2].split(">")[0]
				check = check + 1
	if check != 6:
		sys.exit("Required data not provided for shift. Data needed:\nes_h\nkappa\nes_n_x\nes_N_y\nes_N_z\nnumblobs\n")
	

def read_blob_blocks(lines, blobs, num_blobs):
	for line in lines:
		if line == "<system>\n":
			system_index = lines.index("<system>\n")
			break
	
	check = 0
	i = system_index
	while check < int(myparams.num_blobs):
		i = i + 1
		line = lines[i]
		if line != "\n":
			sline = line.split()
			if sline[0][1:] == "blob>":
				while sline[0] != "</blob>":
					i = i + 1
					line = lines[i]
					if line != "\n":
						sline = line.split()
						if sline[0][1:] == "state":
							blobs[check].state = sline[2].split(">")[0]
						elif sline[0][1:] == "nodes":
							blobs[check].nodes = sline[2].split(">")[0]
						elif sline[0][1:] =="scale":
							blobs[check].scale = sline[2].split(">")[0]
		
				check = check + 1
		
				
if len(sys.argv) != 5:
	sys.exit("Usage: python move_into_box.py [INPUT FFEA FILE] [x] [y] [z]")

#Opening ffea script, getting relevant data and filenames
inputffea = open(sys.argv[1], "r")
lines = inputffea.readlines()
translate = vector3(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))
myparams = params()
read_param_block(lines, myparams)

blobs = [blob() for i in range(int(myparams.num_blobs))]
read_blob_blocks(lines, blobs, myparams.num_blobs)		


#Calculate centroid of entire system of blobs
path = sys.argv[1].rfind("/")
if path == -1:
	path = "./"
else:
	path = sys.argv[1][0:path] + "/"

inputnode = [open(path + blobs[i].nodes, "r") for i in range(int(myparams.num_blobs))]
num_nodes = 0
blobs_centroid = vector3D(0.0,0.0,0.0)
for i in range(int(myparams.num_blobs)):
	lines = inputnode[i].readlines()
	num_surface_nodes = int(lines[2].split()[1])
	num_interior_nodes = int(lines[3].split()[1])
	num_nodes = num_nodes + int(lines[1].split()[1])
	surface_index = lines.index("surface nodes:\n")
	for j in range(surface_index + 1,surface_index+num_surface_nodes):
		this_line = lines[j].split()
		blobs_centroid.x = blobs_centroid.x + float(this_line[0])
		blobs_centroid.y = blobs_centroid.y + float(this_line[1])
		blobs_centroid.z = blobs_centroid.z + float(this_line[2])

	interior_index = lines.index("interior nodes:\n")
	for j in range(interior_index + 1, interior_index+num_interior_nodes):
		this_line = lines[j].split()
		blobs_centroid.x = blobs_centroid.x + float(this_line[0])
		blobs_centroid.y = blobs_centroid.y + float(this_line[1])
		blobs_centroid.z = blobs_centroid.z + float(this_line[2])

blobs_centroid.x = blobs_centroid.x / num_nodes
blobs_centroid.y = blobs_centroid.y / num_nodes
blobs_centroid.z = blobs_centroid.z / num_nodes

#Translation vector
print translate.x, translate.y, translate.z

#Output new node files

for i in range(int(myparams.num_blobs)):
	root_name = blobs[i].nodes[0:-5]
	print root_name
	outputnode = open(root_name + "_translated.node", "w")
	inputnode[i].seek(0)
	lines = inputnode[i].readlines()
	for j in range(5):
		outputnode.write(lines[j])

	num_surface_nodes = int(lines[2].split()[1])
	num_interior_nodes = int(lines[3].split()[1])
	for j in range(5, 5 + num_surface_nodes):
		sline = lines[j].split()
		new_line = str(float(sline[0]) + translate.x) + " " + str(float(sline[1]) + translate.y) + " " + str(float(sline[2]) + translate.z) + "\n"
		outputnode.write(new_line)
	outputnode.write(lines[5 + num_surface_nodes])
	for j in range(5 + num_surface_nodes + 1, 5 + num_surface_nodes + 1 + num_interior_nodes):
		sline = lines[j].split()
		new_line = str(float(sline[0]) + translate.x) + " " + str(float(sline[1]) + translate.y) + " " + str(float(sline[2]) + translate.z) + "\n"
		outputnode.write(new_line)
