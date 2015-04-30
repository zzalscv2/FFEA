import os
import sys
import errno
import string
import datetime
from time import sleep

class param_instruction:

	# For numerical parameters that take a list of values to try out
	def __init__(self, list_of_values):
		self.list_vary = True
		self.list_of_values = list_of_values
		if len(list_of_values) == 1:
			# If there is only one value in the list, then this parameter does not vary
			self.does_vary = False
			self.value = 0
		else:
			# If there are more values in the list, parameter will cycle through each one on each increment
			self.does_vary = True
			self.start_value = 0
			self.end_value = len(list_of_values) - 1
			self.incr = 1
			self.value = self.start_value

	# Only used if this is a variable parameter
	def increment(self):
		if self.value > self.end_value:
			raise "Out of range"

		self.value += self.incr

	def at_end(self):
		if self.value >= self.end_value:
			return True
		else:
			return False

	def reset(self):
		self.value = self.start_value

	# return the current value of this parameter
	def get_value_as_string(self):
		if self.list_vary == True:
			return str(self.list_of_values[self.value])
		else:
			return str(self.value)

	# returns True if this parameter is numerical and varies
	def varies(self):
		return self.does_vary

class walrus_job:
	def __init__(self, walrus_input, arc1_script, index):
		self.walrus_input = walrus_input
		self.arc1_script = arc1_script
		self.index = index

	def submit(self):
		command = "qsub " + self.arc1_script + " > __shit_manager_qsub_output__"
		print command
		os.system(command + "\n")
		qsub_out = open("__shit_manager_qsub_output__", "r")
		line = qsub_out.readline()
		qsub_out.close()
		line_split = line.split()
		self.job_id = line_split[2]
		print "Submitted with job id = " + self.job_id
		self.complete = False

	def status(self):
		os.system("qstat | grep '" + self.job_id + "' > __shit_manager_qstat_output__")
		qstat_out = open("__shit_manager_qstat_output__", "r")
		line = qstat_out.readline()
		qstat_out.close()
		if line == '':
			return 'Not in list'
		line_split = line.split()
		return line_split[4]

	def get_index(self):
		return self.index

	def get_job_id(self):
		return self.job_id

	def is_complete(self):
		return self.complete

	def has_failed(self):
		# check error file for possible causes of failure
		error_fname = "walrus_" + str(self.index) + ".sh.o" + self.job_id
		print "Checking error file (" + error_fname + ") for clues..."
		os.system("cat " + error_fname + " | grep -i 'Error' > __shit_manager_error_check_output__")
		error_check = open("__shit_manager_error_check_output__", "r")
		line = error_check.readline()
		error_check.close()
		if line == '':
			print "No errors..."
		else:
			print "Errors found. Job has failed."
			return True

		os.system("cat " + error_fname + " | grep -i 'segmentation fault' > __shit_manager_error_check_output__")
		error_check = open("__shit_manager_error_check_output__", "r")
		line = error_check.readline()
		error_check.close()
		if line == '':
			print "No segmentation faults..."
		else:
			print "Segmentation fault found. Job has failed."
			return True

		print "No errors or seg faults found. Job " + str(self.index) + " is complete!"
		self.complete = True
		return False

	def resubmit(self):
		# Half the timestep value in the input file
		print "Halving the time step..."
		os.system("sh half_dt.sh " + self.walrus_input + "\n")

		# resubmit this job
		print "Resubmitting..."
		self.submit()

# Creates a new directory with the given name, suppressing the "already exists" error (but allowing all other errors to be raised)
def create_directory(dirname):
    try:
        os.makedirs(dirname)
    except OSError, exception:
        if exception.errno != errno.EEXIST:
            raise


# Instructions to the manager about which vars to vary and through what range
WALRUS_COMMAND = "~/walrus_6thDec2012/bin/walrus"
blob_stem_name = param_instruction(["~/walrus_6thDec2012/proteins/myoglobin"])
Nx = param_instruction([10])
Ny = param_instruction([10])
Nz = param_instruction([10])
es_h = param_instruction([3])
kappa = param_instruction([1e9])
scale = param_instruction([1.0])
margin_x = param_instruction([0.0])
margin_y = param_instruction([5e-9])
margin_z = param_instruction([0.0])
num_blobs = param_instruction([20])
vdw_eps = param_instruction([1e17, 2e18, 4e18, 6e18, 8e18, 1e19, 2e19])
bulk_and_shear_mod = param_instruction(["20.0e6 10.0e6", "166.0e6 10.0e6", "300.0e6 10.0e6", "1000.0e6 10.0e6", "166.0e6 70.0e6", "300.0e6 70.0e6", "1000.0e6 70.0e6", "300.0e6 100.0e6", "1000.0e6 100.0e6", "1000.0e6 300.0e6"])
dt = param_instruction([1e-14])

instructions = [blob_stem_name, Nx, Ny, Nz, es_h, kappa, scale, margin_x, margin_y, margin_z, num_blobs, vdw_eps, bulk_and_shear_mod, dt]

# Get the list of parameters that vary
vary_params = []
for param in instructions:
	if param.varies() == True:
		vary_params.append(param)

# Create the list of all input combinations due to these varying parameters
input_list = []
done = False
while done == False:
	# Create the input string for this combination of parameters and add to list
	input = ''
	for instr in instructions:
		input += ' ' + instr.get_value_as_string()
	input_list.append(input)

	i = 0
	for param in vary_params:
		if param.at_end() == True:
			param.reset()
			i += 1
		else:
			param.increment()
			break
	if i == len(vary_params):
		done = True

print "List of input strings to be used:"
print "\n".join(input_list)
print "\nThis will require", len(input_list), "simulations."

# Create all the necessary directories each containing its associated walrus input file
print "Creating directories, walrus input scripts and arc1 run scripts..."
now = datetime.datetime.now()
current_dir = os.getcwd()
print "current directory = " + current_dir
master_dir = current_dir + "/run_data_folders_" + str(now).replace(" ", "_")
create_directory(master_dir)
input_index = 0
job_list = []
for input in input_list:

	# Get the unique name for this set of parameters
	sim_name = input.replace(" ", "_").replace("..", "dd").replace("/", "fs").replace("~", "t")

	# Make the directory
	folder_name = master_dir + "/run_" + str(input_index) + sim_name
	print "Creating directory " + folder_name
	create_directory(folder_name)

	# Make the walrus input file
	walrus_input_file = str(input_index) + sim_name + ".walrus"
	walrus_full_path = folder_name + "/" + str(input_index) + sim_name + ".walrus"
	print "Creating walrus input file " + walrus_full_path
	print "python super_tile_plus.py " + input + " > " + walrus_full_path
	os.system("python super_tile_plus.py " + input + " > " + walrus_full_path + "\n")

	# Make the arc1 run script
	arc1_script_fname = folder_name + "/" + "walrus_" + str(input_index) + ".sh"
	print "Creating arc1 run script " + arc1_script_fname
	arc1_script = ""
	arc1_script += "#$ -l h_rt=48:00:00\n"
#	arc1_script += "#$ -pe smp 1\n"
#	arc1_script += "#$ -l cputype=amd\n"
#	arc1_script += "#$ -l placement=optimal\n"
	arc1_script += "#$ -cwd -V\n"
	arc1_script += "cd " + folder_name + "\n"
	arc1_script += "export OMP_NUM_THREADS=1\n"
	arc1_script += WALRUS_COMMAND + " " + walrus_input_file + "\n"

	asfile = open(arc1_script_fname, "w")
	asfile.write(arc1_script)
	asfile.close()

	job_list.append(walrus_job(walrus_full_path, arc1_script_fname, input_index))
	input_index += 1
print "...done."

# Now qsub all of these scripts
print "qsubbing all jobs..."
for job in job_list:
	job.submit()
print "...done."

# Now monitor the jobs, checking for and resubmitting failed jobs (with decreased time steps)
print "Monitoring job status:"
jobs_running = True
while jobs_running == True:
	
	# Print job info
	print "Job index | Job ID | Status | Completed"
	for job in job_list:
		print job.get_index(), job.get_job_id(), job.status(), job.is_complete()

	# Check for failed jobs and resubmit them
	for job in job_list:
		if job.is_complete() == False:
			if job.status() == 'Not in list':
				print "Job " + str(job.get_index()) + " is not in the job list, but has not completed. Investigating..."
				if job.has_failed() == True:
					print "Job " + str(job.get_index()) + " has failed."
					job.resubmit()
		sleep(1)
