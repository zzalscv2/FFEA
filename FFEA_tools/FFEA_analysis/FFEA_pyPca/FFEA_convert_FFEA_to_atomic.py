import sys, os
import FFEA_traj

if len(sys.argv) < 8:
	sys.exit("Usage: python " + sys.argv[0] + " -traj [FFEA trajectory (.out)] -scale [FFEA scale (inverts to work in metres)] -out [OUT fname] -format[mdcrd/pdb] -blob [blob_number (optional)] -frames [num_frames (optional)]")

# Get arguments
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0])) + "/"
traj_fname_set = 0
out_fname_set = 0
format_set = 0
scale_set = 0;
blob_num = -1
num_frames = 100000
for i in range(1, len(sys.argv), 2):
	if sys.argv[i].strip() == "-traj":
		traj_fname = sys.argv[i + 1]
		traj_fname_set = 1
		print "traj_fname = " + traj_fname 
	elif sys.argv[i].strip() == "-out":
		out_fname = sys.argv[i + 1]
		out_fname_set = 1
		print "out_fname = " + out_fname
	elif sys.argv[i].strip() == "-format":
		format = sys.argv[i + 1]
		if format != "pdb" and format != "mdcrd":
			sys.exit("Error. Unexpected format '" + format + "'\nAccepted formats:\nmdcrd\npdb\n")
		format_set = 1
		print "format = " + format
	elif sys.argv[i].strip() == "-scale":
		scale = float(sys.argv[i + 1])
		scale_set = 1
		print "scale = " + str(scale)
	elif sys.argv[i].strip() == "-blob":
		blob_num = int(sys.argv[i + 1])
		print "blob_num = " + str(blob_num)
	elif sys.argv[i].strip() == "-frames":
		num_frames = int(sys.argv[i + 1])
		print "num_frames = " + str(num_frames)
	else:
		sys.exit("Error. Unexpected argument '" + sys.argv[i] + "'\n")

if traj_fname_set == 0 or out_fname_set == 0 or scale_set == 0 or format_set == 0:
	sys.exit("Error. '-traj', '-out', 'format' and '-scale' must all be set.\n")

if len(out_fname.split(".")) == 2:
	out_basename = out_fname.split(".")[0]
else:
	out_basename = out_fname

# Make a new trajectory from blob_num (if necessary)
if blob_num != -1:
	final_traj_fname = FFEA_traj.make_single_blob_traj(traj_fname, 0, num_frames)
else:
	final_traj_fname = traj_fname

traj = FFEA_traj.FFEA_traj(final_traj_fname, num_frames, 0, num_frames, 1)

# Extract first frame from trajectory
first_frame_fname = traj.write_frame_to_file(0)

# Convert first frame to a pseudo .pdb file
first_frame_fname_pdb = first_frame_fname.split(".")[0] + ".pdb"
os.system(scriptdir + "../../FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory/FFEA_convert_traj_to_pdb " + first_frame_fname + " " + first_frame_fname_pdb + " " + str(num_frames) + " " + str(scale))

# Convert entire trajectory to a pdb trajectory
if format == "mdcrd":
	final_traj_fname_pdb = final_traj_fname.split(".")[0] + ".mdcrd"
	os.system(scriptdir + "../../FFEA_analysis/PCA/pyPca/FFEA_convert_traj_to_mdcrd " + final_traj_fname + " " + final_traj_fname_pdb + " " + str(num_frames) + " " + str(scale))

elif format == "pdb":
	final_traj_fname_pdb = final_traj_fname.split(".")[0] + ".pdb"
	os.system(scriptdir + "../../FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory/FFEA_convert_traj_to_pdb " + final_traj_fname + " " + final_traj_fname_pdb + " " + str(num_frames) + " " + str(scale))

# Tell the user what they have
print "\nOrignal trajectory - " + traj_fname
if blob_num != -1:
	print "NEW Single blob trajectory - " + final_traj_fname
print "NEW first frame trajectory - " + first_frame_fname
print "NEW first frame pdb (send to pcazip as topology)- " + first_frame_fname_pdb
print "NEW trajectory (send to pcazip as trajectory)- " + final_traj_fname_pdb
print "\n"
