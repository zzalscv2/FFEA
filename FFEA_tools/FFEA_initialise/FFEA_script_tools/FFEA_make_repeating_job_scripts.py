import sys, os
import glob

if len(sys.argv) != 5:
	sys.exit("Usage: python " + sys.argv[0] + " [Directory containing jobs] [job_time (hours)] [num processors per job] [path to ffea binary]")

# Get args
script_dir = os.path.abspath(sys.argv[1])
job_time = int(sys.argv[2])
num_proc = int(sys.argv[3])
ffea_loc = os.path.dirname(os.path.abspath(sys.argv[4])) + "/ffea"

# Get list of ffea scripts
script_fnames = glob.glob(script_dir + "/*job*.ffea")

if len(script_fnames) == 0:
	sys.exit("No jobs to submit (or haven't been through 'FFEA_make_repeating_ffea_scripts.py)")

# Make job submission scripts
for script_fname in script_fnames:
	job_fname = os.path.splitext(script_fname)[0] + ".sh"
	fout = open(job_fname, "w")
	fout.write("#! /bin/bash\n")
	fout.write("#$ -cwd -V\n")
	fout.write("#$ -l h_rt=%d:00:00\n" % (job_time))
	fout.write("#$ -pe smp %d\n\n" % (num_proc))
	fout.write(ffea_loc + " " + script_fname)
	fout.close()
