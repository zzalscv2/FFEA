import sys, os
import numpy as np
import FFEA_node, FFEA_kinetic_map, FFEA_pdb

if len(sys.argv) < 4:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .node fname] [OUTPUT .node/.pdb fname] [INPUT ffea .map fname]\n")

# Get args
innode = sys.argv[1]
outfname = sys.argv[2]
inmap = sys.argv[3]

if len(sys.argv) == 5:
	pdbtop = sys.argv[4]
	pdb = FFEA_pdb.FFEA_pdb(pdbtop)
	
# Get nodes
input_nodes = FFEA_node.FFEA_node(innode)

# Get map
kinetic_map = FFEA_kinetic_map.FFEA_kinetic_map(inmap)

# Apply matrix!
output_nodes = kinetic_map.apply_sparse(input_nodes)

# Print to file
base, ext = os.path.splitext(outfname)


# Finding these values would take to long, so set all to surface as this is just a test script
num_nodes = len(output_nodes)
num_surface_nodes = len(output_nodes)
num_interior_nodes = 0

if ext == ".node":
	fout = open(outfname, "w")
	fout.write("ffea node file\n")
	fout.write("num_nodes %d\n" % (num_nodes))
	fout.write("num_surface_nodes %d\n" % (num_surface_nodes))
	fout.write("num_interior_nodes %d\n" % (num_interior_nodes))

	# Surface nodes
	fout.write("surface nodes:\n")
	for i in range(num_surface_nodes):
		fout.write("%6.3f %6.3f %6.3f\n" % (output_nodes[i][0], output_nodes[i][1], output_nodes[i][2]))

	# Interior nodes
	fout.write("interior nodes:\n")
	for i in range(num_surface_nodes, num_nodes, 1):
		fout.write("%6.3f %6.3f %6.3f\n" % (output_nodes[i][0], output_nodes[i][1], output_nodes[i][2]))

	fout.close()
	
elif ext == ".pdb":
	
	# Edit pdbtop object
	pdb.blob = [pdb.blob[0]]
	pdb.num_blobs = 1
	
	pdb.blob[0].frame = [pdb.blob[0].frame[0]]
	pdb.blob[0].num_frames = 1
	pdb.num_frames = 1
	
	if pdb.blob[0].num_atoms != num_nodes:
		sys.exit("Error. PDB base has %d nodes, we have coordinates for %d nodes." % (pdb.blob[0].num_atoms, num_nodes))
		
	pdb.blob[0].frame[0].pos = output_nodes
	pdb.write_to_file(outfname)