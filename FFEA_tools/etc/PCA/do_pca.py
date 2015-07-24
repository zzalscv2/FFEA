import os, sys

if len(sys.argv) != 3 and len(sys.argv) != 4 and len(sys.argv) != 5:
	sys.exit("Usage: python do_pca.py [INPUT WALRUS TRAJECTORY FILE] [NUM MODES TO OUTPUT] OPTIONAL {[BLOB NUMBER]} {[NAME OF OUTPUT PCZ FILE]}")

trajectory_fname = sys.argv[1]
trajectory_fname_stripped = os.path.splitext(os.path.basename(trajectory_fname))[0]
num_modes = int(sys.argv[2])

blob_number = -1
output_pcz = "__do_pca_temp__.pcz"

if len(sys.argv) >= 4:
	if len(sys.argv[3]) == 1:
		blob_number = sys.argv[3]
		blob_traj_fname = "__blob_traj__.traj"
	else:
		output_pcz = sys.argv[3]


if len(sys.argv) == 5:
	output_pcz = sys.argv[4]


print "trajectory filename = " + trajectory_fname
print "trajectory filename stripped = " + trajectory_fname_stripped
print "num_modes = " + str(num_modes)
if blob_number == -1:
	print "blob_number = No single blob selected"
else:
	print "blob_number = " + blob_number

sys.exit()
#Get a single blob if necessary
if blob_number != -1:
	print "Extracting single blob from file..."
	os.system("python extract_single_blob.py " + trajectory_fname + " " + blob_traj_fname + " " + blob_number) 
	input_filename = blob_traj_fname
else:
	input_filename = trajectory_fname

# get number of nodes from traj file
os.system("head -n4 " + input_filename + " | tail -n1 > __do_pca_temp_num_nodes__\n")

num_nodes_file = open("__do_pca_temp_num_nodes__", "r")
num_nodes = int(num_nodes_file.readline())
num_nodes_file.close()
print "Number of nodes =", num_nodes

os.system("./walrus_to_x " + input_filename + " __do_pca_temp__.x\n")
os.system("python first_frame_to_pdb.py " + input_filename + " __do_pca_temp__.pdb\n")

if os.path.exists(output_pcz):
    os.remove(output_pcz)

os.system("./pcazip-4.1/bin/pcazip -i __do_pca_temp__.x -o " + output_pcz + " -n " + str(num_nodes) + "\n")

for i in range(1, num_modes+1):
	mode_name = trajectory_fname_stripped + "_anim_" + str(i) + ".pdb"

	if os.path.exists(mode_name):
	    os.remove(mode_name)

	print "Making " + mode_name
	os.system("./pcazip-4.1/bin/pczdump -i " + output_pcz + " -anim " + str(i) + " -pdb __do_pca_temp__.pdb -o " + mode_name + "\n")

	walrus_traj_name = trajectory_fname_stripped + "_anim_" + str(i) + "_walrus_format_trajectory.out"
	print "Making " + walrus_traj_name
	os.system("python convert_pdb_to_walrus.py " + mode_name + " " + walrus_traj_name + "\n")
