import os, sys, math

def mat_mult(mat, vec):
	new_vec = [0.0,0.0,0.0]
	for i in range(3):
		for j in range(3):
			new_vec[i] = new_vec[i] + mat[i][j] * vec[j]
	return new_vec

def vol_translate(trans):
	for i in range(num_nodes):
		for j in range(3):
			nodes[i][j] = nodes[i][j] + trans[j]
			

def vol_rotate(rot, cent):
	for j in range(3):
		rot[j] = math.radians(rot[j])

	rot_mat = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
	rot_mat[0][0] = math.cos(rot[1]) * math.cos(rot[2])
	rot_mat[0][1] = math.cos(rot[0]) * math.sin(rot[2]) + math.cos(rot[2]) * math.sin(rot[0]) * math.sin(rot[1])
	rot_mat[0][2] = math.sin(rot[0]) * math.sin(rot[2]) - math.cos(rot[0]) * math.sin(rot[1]) * math.cos(rot[2])
	rot_mat[1][0] = -math.cos(rot[1]) * math.sin(rot[2])
	rot_mat[1][1] = math.cos(rot[0]) * math.cos(rot[2]) - math.sin(rot[0]) * math.sin(rot[1]) * math.sin(rot[2])
	rot_mat[1][2] = math.sin(rot[0]) * math.cos(rot[2]) + math.cos(rot[0]) * math.sin(rot[1]) * math.sin(rot[2])
	rot_mat[2][0] = math.sin(rot[1])
	rot_mat[2][1] = -math.sin(rot[0]) * math.cos(rot[1])
	rot_mat[2][2] = math.cos(rot[0]) * math.cos(rot[1])

	node_vec = [0.0, 0.0, 0.0]
	for i in range(num_nodes):
		for j in range(3):
			node_vec[j] = nodes[i][j] - cent[j]

		new_node_vec = mat_mult(rot_mat, node_vec)

		for j in range(3):
			nodes[i][j] = new_node_vec[j] + cent[j]
			
def write_output():
	for line in lines:
		outputvol.write(line)
		if line == "points\n":
			outputvol.write(str(num_nodes) + "\n")
			break
	for i in range(num_nodes):
		line = str(nodes[i][0]).rjust(13) + " " + str(nodes[i][1]).rjust(13) + " " + str(nodes[i][2]).rjust(13) + "\n"
		outputvol.write(line)



if len(sys.argv) != 3:
	sys.exit("Usage: python transform_vol.py [INPUT VOL FILE] [OUTPUT .VOL FILE]")

inputvol = open(sys.argv[1], "r")
outputvol = open(sys.argv[2], "w")

lines = inputvol.readlines()
i = lines.index("points\n")
num_nodes = int(lines[i+1])
nodes = []
for j in range(i+2, i + 2 + num_nodes, 1):
	sline = lines[j].split()
	s = [float(sline[0]), float(sline[1]), float(sline[2])]
	nodes.append(s)

done = 0
while(done == 0):
	transform = raw_input("Would you like to translate (t), rotate (r) or finish (f)?:")
	if(transform == "t"):
		translation = [0.0,0.0,0.0]
		translation[0] = input("Enter x translation:")
		translation[1] = input("Enter y translation:")
		translation[2] = input("Enter z translation:")
		vol_translate(translation)
	elif(transform == "r"):
		rotation = [0.0,0.0,0.0]
		center = [0.0,0.0,0.0]
		rotation[0] = input("Enter angle (in degrees) for x rotation:")
		rotation[1] = input("Enter angle (in degrees) for y rotation:")
		rotation[2] = input("Enter angle (in degrees) for z rotation:")
		center[0] = input("Enter center of rotation, x:")
		center[1] = input("Enter center of rotation, y:")
		center[2] = input("Enter center of rotation, z:")

		vol_rotate(rotation, center)
	elif(transform == "f"):
		write_output()
		done = 1
	else:
		print "Please enter either 'translate', 'rotate' or 'finish'\n"
