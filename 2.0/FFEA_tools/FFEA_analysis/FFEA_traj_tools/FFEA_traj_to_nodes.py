import sys, os
import FFEA_traj

if len(sys.argv) < 4:
	sys.exit("python " + sys.argv[0] + " [INPUT FFEA trajectory file] [num_frames_to_read] {[Blob index] [Conformation index] [INPUT .node fname] [OUTPUT .node fname] [FFEA scale]} * num_blobs")

# Get args
traj_fname = sys.argv[1]
num_frames_to_read = int(sys.argv[2])

blob_index = []
conformation_index = []
innode_fname = []
outnode_fname = []
scale = []
for i in range(3, len(sys.argv), 5):
	blob_index.append(int(sys.argv[i]))
	conformation_index.append(int(sys.argv[i + 1]))
	innode_fname.append(sys.argv[i + 2])
	outnode_fname.append(sys.argv[i + 3])
	scale.append(float(sys.argv[i + 4]))

num_blobs = len(innode_fname)

# Read traj
traj = FFEA_traj.FFEA_traj(traj_fname, num_frames_to_read, 0, num_frames_to_read, 1)
if traj.total_num_blobs != num_blobs:
	sys.exit("Error. 'total_num_blobs' from trajectory not equal to num file pairs provided.")

# Get nodes and stuff
num_nodes = []
num_surface_nodes = []
num_interior_nodes = []
for fname in innode_fname:
	fin = open(fname, "r")
	fin.readline()
	num_nodes.append(int(fin.readline().split()[1]))
	num_surface_nodes.append(int(fin.readline().split()[1]))
	num_interior_nodes.append(int(fin.readline().split()[1]))
	fin.close()

for i in range(len(outnode_fname)):
	fout = open(outnode_fname[i], "w")
	fout.write("ffea node file\n")
	fout.write("num_nodes %d\n" % (num_nodes[i]))
	fout.write("num_surface_nodes %d\n" % (num_surface_nodes[i]))
	fout.write("num_interior_nodes %d\n" % (num_interior_nodes[i]))

	# Surface nodes
	fout.write("surface nodes:\n")
	for j in range(num_surface_nodes[i]):
		x = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][0] / scale[i]
		y = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][1] / scale[i]
		z = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][2] / scale[i]
		fout.write("%5.2e %5.2e %5.2e\n" % (x, y, z))

	# Interior nodes
	fout.write("interior nodes:\n")
	for j in range(num_surface_nodes[i], num_nodes[i], 1):
		x = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][0] / scale[i]
		y = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][1] / scale[i]
		z = traj.blob[blob_index[i]][conformation_index[i]].frame[-1].node_pos[j][2] / scale[i]
		fout.write("%5.2e %5.2e %5.2e\n" % (x, y, z))

	fout.close()
