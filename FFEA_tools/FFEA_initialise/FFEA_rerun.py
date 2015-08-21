import sys, os, time
import FFEA_script, FFEA_traj
if len(sys.argv) != 4:
	sys.exit("Usage: python FFEA_rerun.py [INPUT .ffea file] [num_frames to delete each restart] [num_frames required]")

# Get args
script_fname = os.path.abspath(sys.argv[1])
delete_frames = int(sys.argv[2])
num_frames_required = int(sys.argv[3])

# Do stuff with args
os.chdir(os.path.dirname(script_fname))
script = FFEA_script.FFEA_script(script_fname)
if num_frames_required > script.params.num_steps / script.params.check:
	print "num_frames_required larger than number set in script. Resetting..."
	num_frames_required = script.params.num_steps / script.params.check

# Run endless simulations like a boss
original_script_fname = script_fname
restarts = 0
num_frames_completed = 0
while(True):

	# If on zeroth restart, make sure script will run indefinitely
	if restarts == 0:
		script.params.num_steps = script.params.check * num_frames_required
		script.write_to_file(script_fname.split(".")[0] + "_newnumframes.ffea")

	# If on first restart, script needs changing
	if restarts == 1:
		script.params.restart = 1
		script_fname = script_fname.split(".")[0] + "_restarted.ffea"
		script.write_to_file(script_fname)
		script = FFEA_script.FFEA_script(script_fname)

	# Run sim if possible
	if num_frames_completed < delete_frames and script.params.restart == 1:
		print "Error. Not enough frames to delete for the restart. Fix your simulation parameters."
		break
	else:
		os.system("ffea -l " + str(delete_frames) + " " + script_fname)

	# Completed runs
	restarts += 1

	# How many frames calculated before crash?
	num_frames_completed = FFEA_traj.get_num_frames(script.params.trajectory_out_fname)

	# Do we need to stop?
	if restarts > num_frames_completed:
		print "Restarting more often than frames are being produced. Time to stop."
		print "num_restarts = " + str(restarts)
		print "num_frames_completed = " + str(num_frames_completed)
		break

	if num_frames_completed >= num_frames_required:
		print "Simulation completed!!!!"
		print "num_restarts = " + str(restarts)
		print "num_restarts / num_frames = " + str(float(restarts) / float(num_frames_completed))
		break

		




