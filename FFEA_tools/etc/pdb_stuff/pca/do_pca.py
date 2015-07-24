import os, sys

if len(sys.argv) != 3 and len(sys.argv) != 4:
	sys.exit("Usage: python do_pca.py [INPUT WALRUS TRAJECTORY FILE] [NUM MODES TO OUTPUT] OPTIONAL {[NAME OF OUTPUT PCZ FILE]}")

trajectory_fname = sys.argv[1]
trajectory_fname_stripped = os.path.splitext(os.path.basename(trajectory_fname))[0]
num_modes = int(sys.argv[2])

output_pcz = "__do_pca_temp__.pcz"
if len(sys.argv) == 4:
	output_pcz = sys.argv[3]

print "trajectory filename = " + trajectory_fname
print "trajectory filename stripped = " + trajectory_fname_stripped
print "num_modes = " + str(num_modes)

# get number of nodes from traj file
os.system("head -n3 " + trajectory_fname + " | tail -n1 > __do_pca_temp_num_nodes__\n")

num_nodes_file = open("__do_pca_temp_num_nodes__", "r")
num_nodes = int(num_nodes_file.readline())
num_nodes_file.close()
print "Number of nodes =", num_nodes

os.system("./walrus_to_x " + trajectory_fname + " __do_pca_temp__.x\n")
os.system("python first_frame_to_pdb.py " + trajectory_fname + " __do_pca_temp__.pdb\n")

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
