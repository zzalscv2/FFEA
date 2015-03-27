import sys, os

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT FFEA trajectory file] [num_modes] [variance %]\nIf you want a single blob extracted, run another script first to get the appropriate trajectory file!\n")

intraj = sys.argv[1]
intraj_base, intraj_ext = os.path.splitext(intraj)

num_modes = int(sys.argv[2])

var_percent = float(sys.argv[3])
if var_percent < 0.0 or var_percent > 100.0:
	sys.exit("Error. variance percent must be between 0 and 100.\n")


# Creating trajectory file from FFEA trajectory
temp_x_fname = "__do_pca_temp__.x"
os.system("./ffea_traj_to_x " + intraj +  " " + temp_x_fname + "\n")
print "Completed ffea_traj_to_x: temp_x_fname = " + temp_x_fname + "\n"

# Extracting the first frame of the FFEA trajectory
first_frame_pdb_fname = "__do_pca_temp__.pdb"
os.system("python ffea_first_frame_to_pdb.py " + intraj + " " + first_frame_pdb_fname + "\n")
print "Completed ffea_first_frame_to_pdb: first_frame_pdb_fname = " + first_frame_pdb_fname + "\n"

# PCAZIP time
output_pcz_fname = intraj_base + ".pcz"
os.system("pyPcazip -i " + temp_x_fname + " -o " + output_pcz_fname + " -q " + var_percent)

