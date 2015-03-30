import sys, os
import subprocess
import threading

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [FFEA base fname] [num_jobs] [starting job_num]\nScripts must have gone through 'FFEA_make_multiple_scripts'\n")

script_basename = sys.argv[1].split(".")[0]
script_ext = sys.argv[1].split(".")[1]
num_jobs = int(sys.argv[2])
start_job_num = int(sys.argv[3])

for i in range(start_job_num, start_job_num + num_jobs):
	new_script_fname = script_basename + "_" + str(i) + "." + script_ext
	arg = "/localhome/py09bh/Software/FFEA/FFEA_working_version/bin/ffea " + str(new_script_fname) + " > a.out &\n"
	threading.Thread(target=os.system, args=(arg,)).start()
