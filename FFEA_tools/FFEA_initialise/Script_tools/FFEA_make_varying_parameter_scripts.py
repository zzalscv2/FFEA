import sys, os

if len(sys.argv) != 7:
	sys.exit("Usage: python " + sys.argv[0] + " [orignal_script] [param_name] [min_value] [max_value] [num_jobs] [interval_type(log,linear)]")

# Get args
orig_script_fname = sys.argv[1]
param_name = sys.argv[2]
min_value = float(sys.argv[3])
max_value = float(sys.argv[4])
num_jobs = int(sys.argv[5])
split_type = sys.argv[6]

# Check param_name
if param_name != "stokes_visc":
	sys.exit("Error. Unrecognised parameter '" + param_name + "'")

# Check num_jobs
if num_jobs < 2:
	sys.exit("Error. 'num_jobs' cannot be less than 2.")

# Calc list of values
if split_type == "linear":
	interval = (max_value - min_value) / (num_jobs - 1)
	values = [min_value + i * interval for i in range(num_jobs)]

elif split_type == "log":
	interval = pow(max_value / min_value, 1.0 / (num_jobs - 1))
	values = [min_value * pow(interval, i) for i in range(num_jobs)]
else:
	sys.exit("Error. '" + split_type + "' is not a recognised 'interval_type'.")

# Get input_script
fin = open(orig_script_fname, "r")
lines = fin.readlines()
fin.close()

# Make scripts!
restart = 0
num_blobs = 1
for i in range(num_jobs):
	
	# Get new_script_fname
	new_script_fname = os.path.splitext(orig_script_fname)[0] + "_" + param_name + "_%5.2e.ffea" % (values[i])

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
			new_rvalue = os.path.splitext(rvalue)[0] + "_" + param_name + "_%5.2e.out" % (values[i])
			new_line = line.replace(rvalue, new_rvalue)
			if restart == 1:
				os.system("cp " + rvalue + " " +  new_rvalue)
				
		elif lvalue.strip() == "measurement_out_fname":
			new_rvalue = os.path.splitext(rvalue)[0] + "_" + param_name + "_%5.2e.out" % (values[i])
			new_line = line.replace(rvalue, new_rvalue)
			if restart == 1:
				
				os.system("cp " + os.path.splitext(rvalue)[0] + "_world.out" + " " +  os.path.splitext(new_rvalue)[0] + "_world.out")
				for j in range(num_blobs):
					os.system("cp " + os.path.splitext(rvalue)[0] + "_blob" + str(j) + ".out" + " " + os.path.splitext(new_rvalue)[0] + "_blob" + str(j) + ".out")

		elif lvalue.strip() == param_name:
			new_rvalue = " %5.2e" % (values[i])
			new_line = line.replace(rvalue, new_rvalue)
		else:
			new_line = line
			
		fout.write(new_line)

	# Close file
	fout.close()
