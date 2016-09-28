import sys, os

if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + " [orignal_script] [num_jobs]")

# Get args
orig_script_fname = sys.argv[1]
num_jobs = int(sys.argv[2])

# Get input_script
fin = open(orig_script_fname, "r")
lines = fin.readlines()
fin.close()

# Make scripts!
restart = 0
num_blobs = 1
for i in range(num_jobs):
	
	# Get new_script_fname
	new_script_fname = os.path.splitext(orig_script_fname)[0] + "_job%d.ffea" % (i)

	# Get num_blobs first!
	for line in lines:
		try:
			lvalue, rvalue = line.strip()[1:-1].split("=")
		except ValueError:
			continue

		if lvalue.strip() == "num_blobs":
			num_blobs = int(rvalue)
			break

	# Open file
	fout = open(new_script_fname, "w")

	# Copy and change lines
	for line in lines:
		try:
			lvalue, rvalue = line.strip()[1:-1].split("=")
		except ValueError:
			fout.write(line)
			continue

		if lvalue.strip() == "restart":
			if rvalue == 1:
				restart = 1
			new_line = line

		elif lvalue.strip() == "trajectory_out_fname":
			new_rvalue = os.path.splitext(rvalue)[0] + "_job%d.out" % (i)
			new_line = line.replace(rvalue, new_rvalue)
			if restart == 1:
				os.system("cp " + rvalue + " " + new_rvalue)

		elif lvalue.strip() == "measurement_out_fname":
			new_rvalue = os.path.splitext(rvalue)[0] + "_job%d.out" % (i)
			new_line = line.replace(rvalue, new_rvalue)
			if restart == 1:
				os.system("cp " + os.path.splitext(rvalue)[0] + "_world.out" + " " + os.path.splitext(new_rvalue)[0] + "_world.out")
				for j in range(num_blobs):
					os.system("cp " + os.path.splitext(rvalue)[0] + "_blob%d.out" % (j) + " " + os.path.splitext(new_rvalue)[0] + "_blob%d.out" % (j))
		else:
			new_line = line
		
		fout.write(new_line)

	# Close file
	fout.close()
