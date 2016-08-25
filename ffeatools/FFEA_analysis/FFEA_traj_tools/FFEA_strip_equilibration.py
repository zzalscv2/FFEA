import sys, os, copy
import FFEA_script, FFEA_trajectory, FFEA_measurement
import numpy as np

# Define as a function so we can return eqat if we need it
def strip_equilibration(fname, limit, split=True):

	# Args and load objects
	script = FFEA_script.FFEA_script(fname)
	eqlimit = float(limit) / 100.0
	meas = script.load_measurement()

	print "FFEA_strip_trajectory.py"
	print "\n\tScript that will calculate when you energies equilibrated (within the given limit) and split your simulation into the two phases"

	# We want to check that each individual blob has equilibrated (doing it this way rather than whole system in case we want to be clever with individual blobs in the future)

	# Get the expected values
	keEXP = []
	seEXP = []

	for i in range(script.params.num_blobs):
		top = script.load_topology(i)
		num_nodes = len(top.get_linear_nodes())
		seEXP.append((3 * num_nodes - 6) / 2.0)		# Units of kbT
		if script.blob[i].solver != "CG_nomass":
			keEXP.append(3.0 * num_nodes / 2.0)
		else:
			keEXP.append(0.0)


	# Calculate running averages, stop when we reach equilibration
	num_frames = len(meas.global_meas["Time"])
	kT = script.params.kT
	keavg = [0.0 for i in range(script.params.num_blobs)]
	seavg = [0.0 for i in range(script.params.num_blobs)]

	eqat = -1
	blobs_to_check = range(script.params.num_blobs)
	for t in range(num_frames):

		# Check if finished	
		if blobs_to_check == []:
			eqat = t
			break

		blobs_to_delete = []
		for i in blobs_to_check:
		
			# New averages
			seavg[i] = (seavg[i] * t + meas.blob_meas[i]["StrainEnergy"][t] / kT) / (t + 1)
			if script.blob[i].solver != "CG_nomass":
				keavg[i] = (keavg[i] * t + meas.blob_meas[i]["KineticEnergy"][t] / kT) / (t + 1)
			else:
				keavg[i] = 0.0

			# Check it's within limits (catching zero div errors)
			if np.fabs(seavg[i] - seEXP[i]) / seEXP[i] < eqlimit:
				if script.blob[i].solver != "CG_nomass":
					if np.fabs(keavg[i] - keEXP[i]) / keEXP[i] < eqlimit:
						blobs_to_delete.append(i)

				else:
					blobs_to_delete.append(i)
	
		# Sort out blobs which have already equilibrated
		for i in blobs_to_delete:
			blobs_to_check.remove(i)

	# Check if never equilibrated
	eqtime = meas.global_meas["Time"][eqat]
	print split
	if eqat == -1:
		sys.exit("\n\tYour system never equilibrated with a %5.2f%% limit on the energy values. Please run for longer." % (eqlimit * 100))
	else:
		print "\n\tYour system equilibrated to within a %5.2f%% limit after %d steps, or %e ns." % (eqlimit * 100, eqat, eqtime * 1e9)
		if not split:
			print ("\nRequested not to split. Returning the number of steps to equilibration.")
			return eqat

	# We need to create an equilibration measurement file, and an equilibration trajectory file that can be viewd in the visualiser.
	# It must be exactly the same format as if FFEA had been run in equilibration mode
	# Remaining trajectory and measurement data must have the time coordinates changed and each blob must be moved back to where it started (add option to traj 'start_from' defaulting to zero)
	new_num_frames = num_frames - eqat
	print "\nSplitting system into equilibration (%d frames) and simulation (%d frames) phases...\n" % (eqat, new_num_frames)

	#
	# Traj
	#
	eqtraj = FFEA_trajectory.FFEA_trajectory(script.params.trajectory_out_fname, num_frames_to_read = eqat)
	traj = FFEA_trajectory.FFEA_trajectory(script.params.trajectory_out_fname, start = eqat)

	# If you want to do something clever with the trajectory, like move everything back to the start, then put it here. Triggered by reset flag.
	# You'll have to redo centroids (easy) and RMSDs (not easy, and slow)


	#
	# Meas (get a new one to write out the equilibration stuff too)
	#
	eqmeas = FFEA_measurement.FFEA_measurement()
	eqmeas.time = meas.time
	eqmeas.date = meas.date
	eqmeas.simtype = "Equilibration"

	for i in range(script.params.num_blobs):
		eqmeas.add_empty_blob()
	
	for key in meas.global_meas:
		if key == "Time":
			eqmeas.global_meas[key] = meas.global_meas[key][0:eqat]
			meas.global_meas[key] = meas.global_meas[key][0:new_num_frames]	
		else:
			try:
				eqmeas.global_meas[key] = meas.global_meas[key][0:eqat]
				meas.global_meas[key] = meas.global_meas[key][eqat:]
			except:
				pass

	total_num_nodes = 0
	for i in range(script.params.num_blobs):
		for key in meas.blob_meas[i]:
			try:
				eqmeas.blob_meas[i][key] = meas.blob_meas[i][key][0:eqat]
				meas.blob_meas[i][key] = meas.blob_meas[i][key][eqat:]
			except:
				pass

		for j in range(script.params.num_blobs):
			for key in meas.blob_meas[i]:
				try:
					eqmeas.interblob_meas[i][j][key] = meas.interblob_meas[i][j][key][0:eqat]
					meas.interblob_meas[i][j][key] = meas.interblob_meas[i][j][key][eqat:]
				except:
					pass


	# Now, write all out to files
	mfname = script.params.measurement_out_fname
	mbase, mext = os.path.splitext(script.params.measurement_out_fname)
	meqfname = mbase + "_equilibration" + mext
	mremfname = mbase + "_remainder" + mext

	eqmeas.write_to_file(meqfname, script = script)
	meas.write_to_file(mremfname, script = script)

	tfname = script.params.trajectory_out_fname
	tbase, text = os.path.splitext(script.params.trajectory_out_fname)
	teqfname = tbase + "_equilibration" + text
	tremfname = tbase + "_remainder" + text

	eqtraj.write_to_file(teqfname)
	traj.write_to_file(tremfname)

if __name__ == '__main__':
	try:
		strip_equilibration(sys.argv[1], sys.argv[2], split=sys.argv[3].lower())
	except(IndexError):
		try:
			strip_equilibration(sys.argv[1], sys.argv[2])
		except:
			print "Error calling strip_equilibration function. Command line arguments 'fname' or 'limit' must be incorrect"
			raise
	except:
		print "Error calling strip_equilibration function. Command line arguments must be incorrect"
		raise
