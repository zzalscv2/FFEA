import sys, os
import FFEA_script, FFEA_springs, FFEA_trajectory
import numpy as np
from subprocess import Popen

if len(sys.argv) != 6:
	sys.exit("Usage: python FFEA_make_structures_overlap.py [INPUT .ffea file] [Blob index 1] [Conformation Index 1] [Blob index 2] [Conformation Index 2]")

# Get args
inffea = sys.argv[1]
blob_index = [int(sys.argv[2]), int(sys.argv[4])]
conf_index = [int(sys.argv[3]), int(sys.argv[5])]

# Build an ffea script object
script = FFEA_script.FFEA_script(inffea)

# Now, build an new one for running
outffea = os.path.dirname(os.path.abspath(inffea)) + "/lol.ffea"
script2 = FFEA_script.FFEA_script("")

# Copy params over
script2.params = script.params
script2.params.num_blobs = 0
script2.params.num_conformations = [1,1] 
script2.params.num_states = [1,1]
script2.params.calc_kinetics = 0
script2.params.calc_noise = 0
#script2.params.calc_vdw = 1
#script2.params.vdw_type = "lennard-jones"
#script2.params.vdw_steric_factor = 1e-2
script2.params.trajectory_out_fname = os.path.dirname(os.path.abspath(inffea)) + "/lol_traj.out"
script2.params.measurement_out_basefname = os.path.dirname(os.path.abspath(inffea)) + "/lol_meas.out"

# Now get blobs
for i in range(2):
	script2.add_empty_blob()
	script2.blob[-1].add_empty_conformation()
	script2.blob[-1].conformation[-1].motion_state = script.blob[blob_index[i]].conformation[conf_index[i]].motion_state
	script2.blob[-1].conformation[-1].nodes = script.blob[blob_index[i]].conformation[conf_index[i]].nodes
	script2.blob[-1].conformation[-1].surface = script.blob[blob_index[i]].conformation[conf_index[i]].surface
	script2.blob[-1].conformation[-1].topology = script.blob[blob_index[i]].conformation[conf_index[i]].topology
	script2.blob[-1].conformation[-1].material = script.blob[blob_index[i]].conformation[conf_index[i]].material
	script2.blob[-1].conformation[-1].stokes = script.blob[blob_index[i]].conformation[conf_index[i]].stokes
	script2.blob[-1].conformation[-1].vdw = script.blob[blob_index[i]].conformation[conf_index[i]].vdw
	script2.blob[-1].conformation[-1].pin = script.blob[blob_index[i]].conformation[conf_index[i]].pin

	script2.blob[-1].solver = script.blob[blob_index[i]].solver
	script2.blob[-1].scale = script.blob[blob_index[i]].scale
	
	# Overlap the centroids
	#script2.blob[-1].centroid = np.array([0.0,0.0,0.0])
	script2.blob[-1].centroid = script.blob[blob_index[i]].centroid
	script2.blob[-1].rotation = script.blob[blob_index[i]].rotation

# And add some springs
spring_fname = os.path.dirname(os.path.abspath(inffea)) + "/lol.spring"
script2.spring = spring_fname

# Sort checkpoint files
script2.params.checkpoint_out = os.path.dirname(os.path.abspath(inffea)) + "/lol_checkout.fcp"
script2.params.checkpoint_in = os.path.dirname(os.path.abspath(inffea)) + "/lol_checkin.fcp"
script2.params.calc_vdw = 0
script2.write_to_file(outffea)	
logfile = os.path.dirname(os.path.abspath(inffea)) + "/lol.log"
fout = open(logfile, "w")

# Add some springs into the spring array until completed!
spring_array = []

run = 0
while True:

	#viewer_process = Popen(["python",os.path.expandvars("$FFEAVIEWER") + "/FFEA_viewer.py", "-s", outffea], stdout=fout)
	run += 1
	print("Minimisation Run %d\n\n" % (run))
	line = raw_input("\tPlease enter:\n\t\tPairs of 2 nodes between which you would like to be added as springs\n\t\tPairs of 2 nodes making up a springs you want to delete\n\t\tReturn to continue\n\t\t'q' to finish: ")
	if line.strip() == "q" or line.strip() == "Q":
		print("Minimisation completed!")
		break
	elif line.strip() == "":
		script2.params.num_steps = script2.params.check * 10 * run
		if run == 1:
			script2.params.restart = 0
		else:
			script2.params.restart = 1

		# Sort checkpoint files
		os.system("mv " + os.path.dirname(os.path.abspath(inffea)) + "/lol_checkout.fcp " + os.path.dirname(os.path.abspath(inffea)) + "/lol_checkin.fcp")
		#script2.params.checkpoint_out = os.path.dirname(os.path.abspath(inffea)) + "/lol_checkout.fcp"
		#script2.params.checkpoint_in = os.path.dirname(os.path.abspath(inffea)) + "/lol_checkin.fcp"
		script2.write_to_file(outffea)
		os.system("ffea " + outffea)
		continue
	else:
		sline = line.split()
		if len(sline) % 2 != 0:

			print("Please enter pairs of integers separated by whitespace!")
			run -= 1
			continue	

		try:
			pairs = [[int(sline[i].strip()), int(sline[i + 1].strip())] for i in range(0, len(sline), 2)]
		except:
			print("Enter integers dammit!")
			run -= 1
			continue

		# Remove if already there, else add
		for pair in pairs:
			if pair in spring_array:
				spring_array.remove(pair)
			else:
				spring_array.append(pair)

		# Add springs to script
		line = raw_input("\n\t\tSpring constant?:")
		if line.strip() != "":
			try:
				k = float(line)
			except:
				print("Enter floating point numbers!")
				run -= 1
				continue
		else:
			k = 1e-2

		line = raw_input("\n\t\tEquilibrium length?:")
		if line.strip() != "":
			try:
				l = float(line)
			except:
				print("Enter floating point numbers!")
				run -= 1
				continue
		else:
			l = 0

		springs = FFEA_springs.FFEA_springs("")
		for node_pair in spring_array:
			spring = FFEA_springs.FFEA_spring()
			spring.set_properties(k,l,[0,1],[0,0], node_pair)
			springs.add_spring(spring)

		springs.write_to_file(spring_fname)

		# Update run time
		script2.params.num_steps = script2.params.check * 10 * run
		if run == 1:
			script2.params.restart = 0
		else:
			script2.params.restart = 1

		# Sort checkpoint files
		os.system("mv " + os.path.dirname(os.path.abspath(inffea)) + "/lol_checkout.fcp " + os.path.dirname(os.path.abspath(inffea)) + "/lol_checkin.fcp")
		script2.write_to_file(outffea)

		# End the viewer
		#viewer_process.kill()

		# Start ffea
		os.system("ffea " + outffea)

# Finished minimisation! Make node files from this trajectory
traj = FFEA_trajectory.FFEA_trajectory(script2.params.trajectory_out_fname)
#traj.blob[0][0].write_frame_as_nodes("blob0_overlap.node", traj.num_frames - 1, 1.0 / script2.blob[0].scale)
#traj.blob[1][0].write_frame_as_nodes("blob1_overlap.node", traj.num_frames - 1, 1.0 / script2.blob[1].scale)

traj.blob[0][0].frame[-1].scale(1.0 / script2.blob[0].scale)
traj.blob[1][0].frame[-1].scale(1.0 / script2.blob[1].scale)
traj.blob[0][0].frame[-1].write_to_file("blob0_overlap.node")
traj.blob[1][0].frame[-1].write_to_file("blob1_overlap.node")

# Remove all weirdo files
os.system("rm lol*")
