import sys, os
import glob

if len(sys.argv) != 2:
	sys.exit("Usage: python " + sys.argv[0] + " [Directory containing jobs]")

script_dir = os.path.dirname(os.path.abspath(sys.argv[1]) + "/" + sys.argv[1])

# Get list of job scripts
job_fnames = glob.glob(script_dir + "/*job*.sh")

if len(job_fnames) == 0:
	sys.exit("No jobs to submit (or haven't been through 'FFEA_make_repeating_job_scripts.py)")

for job_fname in job_fnames:
	os.system("qsub " + job_fname)
