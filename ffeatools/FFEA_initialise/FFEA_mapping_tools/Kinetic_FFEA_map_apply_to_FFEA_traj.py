import sys, os
import numpy as np
import FFEA_trajectory, FFEA_kinetic_map, FFEA_pdb

if len(sys.argv) < 4:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT traj fname (.out)] [OUTPUT traj fname (.out/.pdb)] [INPUT ffea .map fname] [INPUT topology .pdb]\n")

# Get args
intraj = sys.argv[1]
outtraj = sys.argv[2]
inmap = sys.argv[3]
intop = ""

if len(sys.argv) == 5:
	intop = sys.argv[4]

# Test files
base, ext = os.path.splitext(outtraj)
if ext == ".pdb":

	# We need a topology i.e. the atom types we are mapping onto
	if intop == "":
		sys.exit("Error. If mapping to a pdb, we need a pdb topology file as a template")
		
	pdbtop = FFEA_pdb.FFEA_pdb(intop)

	
# Get nodes
traj = FFEA_trajectory.FFEA_trajectory(intraj)

# Get map
kinetic_map = FFEA_kinetic_map.FFEA_kinetic_map(inmap)

# Test against target topology if necessary
if intop != "":
	if kinetic_map.num_rows != pdbtop.blob[0].num_atoms and kinetic_map.num_rows != sum(pdbtop.num_atoms):
		sys.exit("Error. Provided topology has %d atoms. Map expects %d target atoms." % (pdbtop.blob[0].num_atoms, kinetic_map.num_rows))

# Apply matrix to all possible blobs!
output_nodes = []   # output_nodes[blob][frame][pos]
for b in traj.blob:
	if b[0].num_nodes == kinetic_map.num_columns:
	
		print "Applying to blob ", traj.blob.index(b)
		output_nodes.append([])
		count = 0
		for f in b[0].frame:
		
			output_nodes[-1].append(kinetic_map.apply_sparse(f))
			count += 1
			sys.stdout.write("\r\t%d frames made out of %d" % (count, len(b[0].frame)))
			sys.stdout.flush()
			
# Print to file

if ext == ".pdb":

	# We'll have to use the original pdb as a template. Make as many pdb blobs as necessary
	outpdb = FFEA_pdb.FFEA_pdb("")
	
	for bnodes in output_nodes:
	
		# Add new blob
		outpdb.add_empty_blob()
		
		# Populate with atoms
		for blob in pdbtop.blob:
			outpdb.blob[-1].atom.append(blob.atom)
		outpdb.blob[-1].num_atoms = len(outpdb.blob[-1].atom)

		# Now add frames
		for f in bnodes:
			outpdb.blob[-1].frame.append(FFEA_pdb.FFEA_pdb_frame())
			outpdb.blob[-1].frame[-1].pos = f * 1e10
			#print outpdb.blob[-1].frame[-1].pos
			outpdb.blob[-1].num_frames += 1
	
	outpdb.num_frames = outpdb.blob[-1].num_frames
	outpdb.write_to_file(outtraj)
			
elif ext == ".out":
		
	print "Currently unavailable. Sorry about that. Fix me if you're a developer!"
