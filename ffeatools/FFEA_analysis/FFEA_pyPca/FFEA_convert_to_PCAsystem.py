# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

import sys, os
import FFEA_trajectory, FFEA_pdb

if len(sys.argv) < 6:
	sys.exit("Usage: python " + sys.argv[0] + " -traj [FFEA trajectory (.ftj)] -out [Output fname] -format[mdcrd/pdb] -blob [blob_number (optional)] -frames [num_frames (optional)]")

# Get arguments
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0])) + "/"
traj_fname_set = 0
out_fname_set = 0
format_set = 0
blob_num = -1
num_frames = 100000
for i in range(1, len(sys.argv), 2):
	if sys.argv[i].strip() == "-traj":
		traj_fname = os.path.abspath(sys.argv[i + 1])
		traj_fname_set = 1
		print("traj_fname = " + traj_fname)
	elif sys.argv[i].strip() == "-out":
		out_fname = os.path.abspath(sys.argv[i + 1])
		out_fname_set = 1
		print("out_fname = " + out_fname)
	elif sys.argv[i].strip() == "-format":
		format = sys.argv[i + 1]
		if format != "pdb" and format != "mdcrd":
			sys.exit("Error. Unexpected format '" + format + "'\nAccepted formats:\nmdcrd\npdb\n")
		format_set = 1
		print("format = " + format)
	elif sys.argv[i].strip() == "-blob":
		blob_num = int(sys.argv[i + 1])
		print("blob_num = " + str(blob_num))
	elif sys.argv[i].strip() == "-frames":
		num_frames = int(sys.argv[i + 1])
		print("num_frames = " + str(num_frames))
	else:
		sys.exit("Error. Unexpected argument '" + sys.argv[i] + "'\n")


scale = 1e-10
if traj_fname_set == 0 or out_fname_set == 0 or format_set == 0:
	sys.exit("Error. '-traj', '-out' and 'format' must all be set.\n")

if len(out_fname.split(".")) == 2:
	out_basename = os.path.splitext(out_fname)[0]
else:
	out_basename = out_fname

final_traj_fname = traj_fname
traj = FFEA_trajectory.FFEA_trajectory(final_traj_fname, num_frames_to_read = num_frames)

# Remove all trajectory except for blob in question, if required
if blob_num != -1:
	traj.set_single_blob(blob_num)

# Extract first frame from trajectory
first_frame_fname = out_basename + "_frame0.ftj"
traj.write_to_file(first_frame_fname, frames=[0,1])

# Convert first frame to a pseudo .pdb file
first_frame_fname_pdb = first_frame_fname.split(".")[0] + ".pdb"

os.system("python " + scriptdir + "../../FFEA_analysis/FFEA_traj_tools/FFEA_convert_traj_to_pdb.py " + first_frame_fname + " " + first_frame_fname_pdb + " " + str(num_frames) + " " + str(scale))

# Convert entire trajectory to a pdb trajectory
if format == "mdcrd":
	final_traj_fname_pdb = out_basename + ".mdcrd"
	os.system(scriptdir + "../../FFEA_analysis/FFEA_pyPca/FFEA_convert_traj_to_mdcrd " + final_traj_fname + " " + final_traj_fname_pdb + " " + str(num_frames) + " " + str(scale))
elif format == "pdb":
	final_traj_fname_pdb = out_basename + ".pdb"
	pdb = FFEA_pdb.FFEA_pdb("")
	pdb.build_from_traj(traj, scale = 1.0 / scale)
	pdb.write_to_file(final_traj_fname_pdb)

# Tell the user what they have
print("\nOrignal trajectory - " + traj_fname)
if blob_num != -1:
	print("NEW Single blob trajectory - " + final_traj_fname)
print("NEW first frame trajectory - " + first_frame_fname)
print("NEW first frame pdb (send to pcazip as topology)- " + first_frame_fname_pdb)
print("NEW trajectory (send to pcazip as trajectory)- " + final_traj_fname_pdb)
print("\n")
