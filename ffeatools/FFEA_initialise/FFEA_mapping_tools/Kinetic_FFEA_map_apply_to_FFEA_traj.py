import sys, os
import numpy as np
import FFEA_trajectory, FFEA_kinetic_map, FFEA_pdb

if len(sys.argv) < 4:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT traj fname (.ftj)] [OUTPUT traj fname (.ftj/.pdb)] [INPUT ffea .map fname] [INPUT topology .pdb]\n")

# Get args
intraj = sys.argv[1]
outtraj = sys.argv[2]
inmap = sys.argv[3]
intop = ""

if len(sys.argv) > 4:
	intop = sys.argv[4]

if len(sys.argv) > 5:
	num_frames_to_read = int(sys.argv[5])

# Test files
base, ext = os.path.splitext(outtraj)
if ext == ".pdb":

	# We need a topology i.e. the atom types we are mapping onto
	if intop == "":
		sys.exit("Error. If mapping to a pdb, we need a pdb topology file as a template")
		
	pdbtop = FFEA_pdb.FFEA_pdb(intop)
	
# Get nodes
traj = FFEA_trajectory.FFEA_trajectory(intraj, num_frames_to_read = num_frames_to_read)

# Get map
kinetic_map = FFEA_kinetic_map.FFEA_kinetic_map(inmap)

# Test against target topology if necessary
if intop != "":
	if kinetic_map.num_rows != sum(pdbtop.num_atoms):
		sys.exit("Error. Provided topology has %d atoms. Map expects %d target atoms." % (pdbtop.blob[0].num_atoms, kinetic_map.num_rows))

	#if kinetic_map.num_rows != pdbtop.blob[0].num_atoms and kinetic_map.num_rows != sum(pdbtop.num_atoms):
	#	sys.exit("Error. Provided topology has %d atoms. Map expects %d target atoms." % (pdbtop.blob[0].num_atoms, kinetic_map.num_rows))

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

#if ext == ".pdb":
#
#	# We'll have to use the original pdb as a template. Make as many pdb blobs as necessary
#	outpdb = FFEA_pdb.FFEA_pdb("")
#	
#	bindex = -1
#	for bnodes in output_nodes:
#	
#		bindex += 1
#
#		# Add new blob
#		outpdb.add_empty_blob()
#		
#		# Populate with atoms
#		#for blob in pdbtop.blob:
#			#outpdb.blob[-1].atom = outpdb.blob[-1].atom + blob.atom #smush any sub-regions of PDBs together
#		outpdb.blob[bindex].atom = pdbtop.blob[bindex].atom
#		#print blob.atom.chain, outpdb.blob[-1].atom.chain
#
#		outpdb.blob[bindex].num_atoms = len(outpdb.blob[bindex].atom)
#
#		# Now add frames
#		for f in bnodes:
#			outpdb.blob[bindex].frame.append(FFEA_pdb.FFEA_pdb_frame())
#			outpdb.blob[bindex].frame[-1].pos = f * 1e10
#			#print outpdb.blob[-1].frame[-1].pos
#			outpdb.blob[bindex].num_frames += 1
#	
#	outpdb.num_frames = outpdb.blob[-1].num_frames
#	outpdb.write_to_file(outtraj)

if ext == ".pdb":
	
	# We'll have to use the original pdb as a template. Make as many pdb blobs as necessary
	outpdb = pdbtop
	outpdb.clear_position_data()

	for bnodes in output_nodes:

		findex = -1
		for f in bnodes:
			findex += 1
			end = 0
			for i in range(outpdb.num_chains):
				start = end
				end = start + outpdb.num_atoms[i]
				outpdb.chain[i].add_empty_frame()
				outpdb.chain[i].frame[findex].pos = f[start:end] * 1e10
	outpdb.num_frames = outpdb.chain[i].num_frames
	outpdb.write_to_file(outtraj)

	# Now simply change the positions using the new mapped ones
	#bindex = -1
	#for bnodes in output_nodes:
	#	bindex += 1
	#	for f in bnodes:

	#		end = 0
	#		for pdbb in outpdb.blob:
	#			start = end
	#			end = start + pdbb.num_atoms
	#			pdbb.frame.append(FFEA_pdb.FFEA_pdb_frame())
	#			pdbb.frame[-1].pos = f[start:end] * 1e10
	
	#outpdb.num_frames = len(bnodes[0])
	#outpdb.write_to_file(outtraj)
		
elif ext == ".out":
		
	print "Currently unavailable. Sorry about that. Fix me if you're a developer!"
