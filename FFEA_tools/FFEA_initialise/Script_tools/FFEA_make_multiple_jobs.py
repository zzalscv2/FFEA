import sys, os

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [FFEA_script (.ffea)] [num_jobs] [start job_num]")

script_fname = sys.argv[1]
num_jobs = int(sys.argv[2])
start_job = int(sys.argv[3])

# Get original script
fin = open(script_fname, "r")
orig_lines = fin.readlines()
fin.close()

# Make new script by altering trajectory/mesurement and script names
for i in range(start_job, start_job + num_jobs):
	new_script_fname = script_fname.split(".")[0] + "_" + str(i) + "." + script_fname.split(".")[1]
	fout = open(new_script_fname, "w")
	for line in orig_lines:
		if "trajectory_out_fname" in line or "measurement_out_fname" in line or "stress_out_fname" in line:
			lvalue, rvalue = line.rstrip()[0:-1].split("=")
			if len(rvalue.rsplit(".", 2)) == 2:
				new_rvalue = rvalue.rsplit(".", 2)[0] + "_" + str(i) + "." + rvalue.rsplit(".", 2)[1]
			else:
				new_rvalue = rvalue + "_" + str(i)

			new_line = lvalue + " = " + new_rvalue + ">" + "\n"
			fout.write(new_line)
		else:
			fout.write(line)
