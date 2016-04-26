import sys, os
import FFEA_trajectory, FFEA_measurement, FFEA_topology, FFEA_pin

def get_path_from_script(path, scriptdir):
	if os.path.isabs(path):
		return path
	else:
		return os.path.abspath(os.path.join(scriptdir, path))

class FFEA_script:

	def __init__(self, fname):

		self.reset()

		# Start reading to test
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found. Returning empty object..."
			self.reset()
			return

		if "<param>\n" not in fin.readlines():
			print "Error. File " + fname  + " not an FFEA script file."
			fin.close()
			return

		fin.close()

		try:
			self.params = self.read_params_from_script(fname)
		except:
			print "Error. Couldn't load <params>...</params>"
			self.reset()
			return

		for i in range(self.params.num_blobs):
			try:
				self.blob.append(self.read_blob_from_file(fname, i, self.params.num_conformations[i]))
			except:
				print "Error. Couldn't load <blob>...</blob> " + str(i)
				self.reset()
				return

		# Get springs
		try:
			self.read_springs_from_file(fname)
		except:
			print "Error. Failed to load <spring>...</spring> "
			self.reset()
			return

	def reset(self):
		self.params = None
		self.blob = []
		self.spring = ""

	def add_blob(self, blob):
		self.blob.append(blob)
		self.params.num_blobs += 1
		self.params.num_conformations.append(1)

	def read_params_from_script(self, fname):

		# Get scriptdir
		scriptdir = os.path.dirname(fname)

		# Open script file
		try:
			fin = open(fname, "r")
		except(IOError):
			print "Error. File " + fname  + " not found."
			return
			
		param_lines = fin.readlines()

		# Close file
		fin.close()

		# Get params block
		param_lines = extract_block_from_lines('param', 0, param_lines)

		# Get some parameters
		params = FFEA_script_params()

		# Assign the params
		for param in param_lines:
			try:
				line = param.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].rstrip()
				rvalue = line.split("=")[1].lstrip()
				
			except:
				print "Error. Could not parse param '" + param + "'"
				return

			params.assign_param(lvalue, rvalue, scriptdir = scriptdir)
		
		# Sort measurement names
		base, ext = os.path.splitext(params.measurement_out_basefname)
		for i in range(params.num_blobs):
			params.measurement_out_fname.append(base + "_blob" + str(i) + ext)
		params.measurement_out_fname.append(base + "_world" + ext)

		# Now, if params have not been correctly initialised, use the measurement world file to get them
		if not params.completed():
			try:
				fin = open(params.measurement_out_fname[-1], "r")

			except(IOError):
				print "Error. File " + params.measurement_out_fname[-1]  + " not found."
				return
			
			
			line = fin.readline().strip()
			while line != "Parameters:":
				if line == "":
					line = fin.readline().strip()
					continue
				elif line.split()[0].strip() == "#":
					print "Old measurement file. Parameters not available here."
					return
				
				else:
					line = fin.readline().strip()

			line = fin.readline().strip()
			while line != "":

				try:
					line = param.strip().replace("<", "").replace(">", "")
					lvalue = line.split("=")[0].rstrip()
					rvalue = line.split("=")[1].lstrip()
				except:
					print "Error. Could not parse param '" + param + "'"
					return

				params.assign_param(lvalue, rvalue, scriptdir = scriptdir)
				line = fin.readline().strip()

			fin.close()
			
		return params


	def read_blob_from_file(self, fname, index, num_conformations):


		# Get scriptdir
		scriptdir = os.path.dirname(fname)

		# Open script file
		try:
			fin = open(fname, "r")
		except(IOError):
			print "Error. File " + fname  + " not found."
			return

		blob_lines = fin.readlines()

		# Close file
		fin.close()

		# Get relevent blob block
		blob_lines = extract_block_from_lines('blob', index, blob_lines)

		# Get some parameters
		blob = FFEA_script_blob()

		# Get conformations first
		for i in range(num_conformations):

			# Get a conformation
			conformation = FFEA_script_conformation()
			
			conformation_lines = extract_block_from_lines('conformation', i, blob_lines)
			for line in conformation_lines:
				try:
					line = line.strip().replace("<", "").replace(">", "")
					lvalue = line.split("=")[0].rstrip()
					rvalue = line.split("=")[1].lstrip()
				except:
					print "Error. Could not parse conformation tag '" + line + "'"
					return None
		
				try:
					if lvalue == "motion_state":
						conformation.motion_state = rvalue
					elif lvalue == "nodes":
						conformation.nodes = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "topology":
						conformation.topology = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "surface":
						conformation.surface = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "material":
						conformation.material = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "stokes":
						conformation.stokes = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "vdw":
						conformation.vdw = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "pin":
						conformation.pin = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "binding_sites":
						conformation.bsites = get_path_from_script(rvalue, scriptdir)
					else:
						print "Unrecognised conformation tag '" + line + "'. Ignoring..."
						continue

				except(IndexError, ValueError):
					print "Error. Couldn't parse conformation tag '" + line + "'"
					return None
			
			blob.conformation.append(conformation)
			blob.num_conformations += 1

		# Now get kinetic stuff (if it's there)
		kinetic_lines = extract_block_from_lines('kinetics', 0, blob_lines)
		
		# States and rates
		for line in kinetic_lines:
			try:
				line = line.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].rstrip()
				rvalue = line.split("=")[1].lstrip()
			except(IndexError):
				continue
			except:
				print "Error. Could not parse blob tag '" + line + "'"
				return None

			try:
				if lvalue == "states":
					blob.states = rvalue
				elif lvalue == "rates":
					blob.rates = rvalue
				else:
					continue

			except(IndexError, ValueError):
				print "Error. Couldn't parse blob tag '" + line + "'"
				return None

		# Now maps
		map_lines = extract_block_from_lines('maps', 0, kinetic_lines)
		blob.map_indices = []
		blob.map = []
		for line in map_lines:
			
			try:
				line = line.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].rstrip()
				rvalue = line.split("=")[1].lstrip()

			except:
				print "Error. Could not parse blob tag '" + line + "'"
				return None

			try:
				blob.map_indices.append([int(lvalue.strip().split("(")[1].replace(")","").split(",")[i]) for i in range(2)])
				blob.map_indices
				blob.map.append(rvalue)

			except:
				print "Error. Could not parse blob tag '" + line + "'"
				return None

		# Now blob general stuff
		for line in blob_lines:
			try:
				line = line.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].rstrip()
				rvalue = line.split("=")[1].lstrip()
			except(IndexError):
				continue
			except:
				print "Error. Could not parse blob tag '" + line + "'"
				return None

			try:
				if lvalue == "solver":
					blob.solver = rvalue
				elif lvalue == "scale":
					blob.scale = float(rvalue)
				elif lvalue == "centroid" or lvalue == "centroid_pos":
					blob.centroid = [float(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
				elif lvalue == "rotation":
					blob.rotation = [float(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
				else:
					continue

			except(IndexError, ValueError):
				print "Error. Couldn't parse blob tag '" + line + "'"
				return None

		return blob

	def read_springs_from_file(self, fname):

		fin = open(fname, "r")
		lines = fin.readlines()
		fin.close()
		
		spring_lines = extract_block_from_lines("springs", 0, lines)
		
		if len(spring_lines) == 0:
			return

		if len(spring_lines) != 1:
			print "Error. Expected only one filename."
			return

		line = spring_lines[0]

		try:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].rstrip()
			rvalue = line.split("=")[1].lstrip()

		except(IndexError, ValueError):
			print "Error. Couldn't parse spring tag '" + line + "'"
			return

		if lvalue == "springs_fname":
			self.spring = rvalue
			
		return

	def write_to_file(self, fname):

		fout = open(fname, "w")
		self.params.write_to_file(fout, fname)
		fout.write("<system>\n")
		for blob in self.blob:
			blob.write_to_file(fout, fname, self.params.calc_kinetics)

		if self.spring != "":
			fout.write("\t<interactions>\n\t\t<springs>\n")
			fout.write("\t\t\t<spring_fname = %s>\n" % (os.path.relpath(self.spring, os.path.dirname(os.path.abspath(fname)))))
			fout.write("\t\t</springs>\n\t</interactions>\n")
		fout.write("</system>")
		fout.close()

	def load_topology(self, blob_index, conformation_index):

		return FFEA_topology.FFEA_topology(self.blob[blob_index].conformation[conformation_index].topology)

	def load_pin(self, blob_index, conformation_index):
	
		return FFEA_pin.FFEA_pin(self.blob[blob_index].conformation[conformation_index].pin)
	
	def load_trajectory(self, frames_to_read = float("inf")):
		
		return FFEA_trajectory.FFEA_trajectory(self.params.trajectory_out_fname, num_frames_to_read = frames_to_read)

	def load_measurement(self, frames_to_read = float("inf")):

		return FFEA_measurement.FFEA_measurement(self.params.measurement_out_basefname, self.params.num_blobs, num_frames_to_read = frames_to_read)

	def add_empty_blob(self):

		self.blob.append(FFEA_script_blob())
		self.params.num_blobs += 1

class FFEA_script_params():
	
	def __init__(self):
		self.restart = 0
		self.dt = 0.0
		self.kT = 0.0
		self.check = 0
		self.num_steps = 0
		self.trajectory_out_fname = ""
		self.measurement_out_basefname = ""
		self.measurement_out_fname = []
		self.vdw_forcefield_params = ""
		self.kinetics_out_fname = ""
		self.binding_site_params = ""
		self.epsilon = 0.0
		self.max_iterations_cg = 0
		self.kappa = 0.0
		self.epsilon_0 = 0.0
		self.dielec_ext = 0.0
		self.calc_stokes = 0
		self.stokes_visc = 0.0
		self.calc_vdw = 0
		self.calc_noise = 0
		self.calc_es = 0
		self.calc_kinetics = 0
		self.kinetics_update = 0
		self.es_update = 0
		self.es_N_x = 0
		self.es_N_y = 0
		self.es_N_z = 0
		self.vdw_type = "steric"
		self.vdw_steric_factor = 0.0
		self.move_into_box = 1
		self.sticky_wall_xz = 0
		self.wall_x_1 = ""
		self.wall_x_2 = ""
		self.wall_y_1 = ""
		self.wall_y_2 = ""
		self.wall_z_1 = ""
		self.wall_z_2 = ""
		self.es_h = 0
		self.num_blobs = 0
		self.num_conformations = []
		self.num_states = []

	def assign_param(self, lvalue, rvalue, scriptdir = "."):

		if lvalue == "restart":
			self.restart = int(rvalue)
		elif lvalue == "dt":
			self.dt = float(rvalue)
		elif lvalue == "rng_seed":
			pass
		elif lvalue == "kT":
			self.kT = float(rvalue)
		elif lvalue == "check":
			self.check = int(rvalue)
		elif lvalue == "num_steps":
			self.num_steps = int(float(rvalue))
		elif lvalue == "trajectory_out_fname":
			self.trajectory_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "measurement_out_fname":
			self.measurement_out_basefname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "vdw_forcefield_params":
			self.vdw_forcefield_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "kinetics_out_fname":
			self.kinetics_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "binding_site_params":
			self.binding_site_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "epsilon":
			self.epsilon = float(rvalue)
		elif lvalue == "max_iterations_cg":
			self.max_iterations_cg = int(rvalue)
		elif lvalue == "kappa":
			self.kappa = float(rvalue)
		elif lvalue == "epsilon_0":
			self.epsilon_0 = float(rvalue)
		elif lvalue == "dielec_ext":
			self.dielec_ext = float(rvalue)
		elif lvalue == "calc_stokes":
			self.calc_stokes = int(rvalue)
		elif lvalue == "calc_kinetics":
			self.calc_kinetics = int(rvalue)
		elif lvalue == "kinetics_update":
			self.kinetics_update = int(rvalue)
		elif lvalue == "stokes_visc":
			self.stokes_visc = float(rvalue)
		elif lvalue == "calc_vdw":
			self.calc_vdw = int(rvalue)
		elif lvalue == "vdw_type":
			self.vdw_type = rvalue
		elif lvalue == "vdw_steric_factor":
			self.vdw_steric_factor = float(rvalue)
		elif lvalue == "calc_noise":
			self.calc_noise = int(rvalue)
		elif lvalue == "calc_es":
			self.calc_es = int(rvalue)
		elif lvalue == "es_update":
			self.es_update = int(rvalue)
		elif lvalue == "es_N_x":
			self.es_N_x = int(rvalue)
		elif lvalue == "es_N_y":
			self.es_N_y = int(rvalue)
		elif lvalue == "es_N_z":
			self.es_N_z = int(rvalue)
		elif lvalue == "move_into_box":
			self.move_into_box = int(rvalue)
		elif lvalue == "sticky_wall_xz":
			self.sticky_wall_xz = int(rvalue)
		elif lvalue == "wall_x_1":
			self.wall_x_1 = rvalue
		elif lvalue == "wall_x_2":
			self.wall_x_2 = rvalue
		elif lvalue == "wall_y_1":
			self.wall_y_1 = rvalue
		elif lvalue == "wall_y_2":
			self.wall_y_2 = rvalue
		elif lvalue == "wall_z_1":
			self.wall_z_1 = rvalue
		elif lvalue == "wall_z_2":
			self.wall_z_2 = rvalue
		elif lvalue == "es_h":
			self.es_h = int(rvalue)
		elif lvalue == "num_blobs":
			self.num_blobs = int(rvalue)
		elif lvalue == "num_conformations":
			self.num_conformations = [int(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
		elif lvalue == "num_states":
			self.num_states = [int(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
		else:
			print "Unrecognised parameter '" + param + "'. Ignoring..."

	# This function tests whether or not there are enough params to form an ffea system
	def completed(self):
		return True
	
	def test_traj(self):

		# Does the trajectory exist and is it usable?
		if not os.path.exists(self.trajectory_out_fname):
			return False
		
		if not os.path.isfile(self.trajectory_out_fname):
			return False
		
		fin = open(self.trajectory_out_fname)
		line = fin.readline().strip()
		fin.close()
		if line != "FFEA_trajectory_file":
			return False
		else:
			return True

	def write_to_file(self, fout, fname):
		
		num_conformations_string = ""
		num_states_string = ""
		for i in range(self.num_blobs):
			num_conformations_string += str(self.num_conformations[i]) + ","
			#num_states_string += str(self.num_states[i]) + ","

		num_conformations_string = num_conformations_string[0:-1]
		#num_states_string = num_states_string[0:-1]
		

		astr = ""
		astr += "<param>\n"
		astr += "\t<restart = %d>\n" % (self.restart)
		astr += "\t<dt = %5.2e>\n" % (self.dt)
		astr += "\t<kT = %5.2e>\n" % (self.kT)
		astr += "\t<check = %d>\n" % (self.check)
		astr += "\t<num_steps = %1.0e>\n" % (self.num_steps)
		astr += "\t<rng_seed = time>\n"
		astr += "\t<trajectory_out_fname = %s>\n" % (os.path.relpath(self.trajectory_out_fname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<measurement_out_fname = %s>\n" % (os.path.relpath(self.measurement_out_basefname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<vdw_forcefield_params = %s>\n" % (os.path.relpath(self.vdw_forcefield_params, os.path.dirname(os.path.abspath(fname))))
		if self.kinetics_out_fname != "":
			astr += "\t<kinetics_out_fname = %s>\n" % (os.path.relpath(self.kinetics_out_fname, os.path.dirname(os.path.abspath(fname))))

		if self.binding_site_params != "":
			astr += "\t<binding_site_params = %s>\n" % (os.path.relpath(self.binding_site_params, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<epsilon = %5.2e>\n" % (self.epsilon)
		astr += "\t<max_iterations_cg = %d>\n" % (self.max_iterations_cg)
		astr += "\t<kappa = %5.2e>\n" % (self.kappa)
		astr += "\t<epsilon_0 = %5.2e>\n" % (self.epsilon_0)
		astr += "\t<dielec_ext = %5.2e>\n" % (self.dielec_ext)
		astr += "\t<calc_stokes = %d>\n" % (self.calc_stokes)
		astr += "\t<stokes_visc = %5.2e>\n" % (self.stokes_visc)
		astr += "\t<calc_vdw = %d>\n" % (self.calc_vdw)
		if self.calc_vdw == 1:
			astr += "\t<vdw_type = %s>\n" % (self.vdw_type)
			if self.vdw_type == "steric" or self.vdw_type == "ljsteric":
				astr += "\t<vdw_steric_factor = %f>\n" % (self.vdw_steric_factor)

		astr += "\t<calc_noise = %d>\n" % (self.calc_noise)
		astr += "\t<calc_es = %d>\n" % (self.calc_es)
		astr += "\t<es_update = %d>\n" % (self.es_update)
		astr += "\t<es_N_x = %d>\n" % (self.es_N_x)
		astr += "\t<es_N_y = %d>\n" % (self.es_N_y)
		astr += "\t<es_N_z = %d>\n" % (self.es_N_z)
		astr += "\t<move_into_box = %d>\n" % (self.move_into_box)
		astr += "\t<sticky_wall_xz = %d>\n" % (self.sticky_wall_xz)
		astr += "\t<wall_x_1 = %s>\n" % (self.wall_x_1)
		astr += "\t<wall_x_2 = %s>\n" % (self.wall_x_2)
		astr += "\t<wall_y_1 = %s>\n" % (self.wall_y_1)
		astr += "\t<wall_y_2 = %s>\n" % (self.wall_y_2)
		astr += "\t<wall_z_1 = %s>\n" % (self.wall_z_1)
		astr += "\t<wall_z_2 = %s>\n" % (self.wall_z_2)
		astr += "\t<es_h = %d>\n" % (self.es_h)
		astr += "\t<num_blobs = %d>\n" % (self.num_blobs)
		astr += "\t<num_conformations = (%s)>\n" % (num_conformations_string)
		#astr += "\t<num_states = (%s)>\n" % (num_states_string)
		astr += "</param>\n\n"	
		fout.write(astr)

class FFEA_script_blob:

	def __init__(self):
		
		self.num_conformations = 0
		self.conformation = []
		self.states = ""
		self.rates = ""
		self.map = []
		self.map_indices = []
		self.solver = ""
		self.scale = 0.0
		self.centroid = None
		self.rotation = None

	def write_to_file(self, fout, fname, calc_kinetics):

		fout.write("\t<blob>\n")
		need_solver = 0;
		for conformation in self.conformation:
			conformation.write_to_file(fout, fname, calc_kinetics)
			if conformation.motion_state == "DYNAMIC":
				need_solver = 1
		
		if calc_kinetics == 1 and self.states != "":
			fout.write("\t\t<kinetics>\n")
			if self.num_conformations != 0:
				fout.write("\t\t\t<maps>\n")
				for i in range(len(self.map)):
					fout.write("\t\t\t\t<map (%d,%d) = %s>\n" % (self.map_indices[i][0], self.map_indices[i][1], os.path.relpath(self.map[i], os.path.dirname(os.path.abspath(fname)))))
				fout.write("\t\t\t</maps>\n")
			
			fout.write("\t\t\t<states = %s>\n" % (os.path.relpath(self.states, os.path.dirname(os.path.abspath(fname)))))
			fout.write("\t\t\t<rates = %s>\n" % (os.path.relpath(self.rates, os.path.dirname(os.path.abspath(fname)))))
			fout.write("\t\t</kinetics>\n")

		if need_solver == 1:
			fout.write("\t\t<solver = %s>\n" % (self.solver))

		fout.write("\t\t<scale = %5.2e>\n" % (self.scale))
		if self.rotation != None:
			fout.write("\t\t<rotation = (")
			for i in range(len(self.rotation)):
				fout.write("%5.2f" % (self.rotation[i]))
				if i != len(self.rotation) - 1:
					fout.write(",")

			fout.write(")>\n")
		if self.centroid != None:
			fout.write("\t\t<centroid = (")
			for i in range(len(self.centroid)):
				fout.write("%5.2f" % (self.centroid[i]))
				if i != len(self.centroid) - 1:
					fout.write(",")

			fout.write(")>\n")
		fout.write("\t</blob>\n")

	def add_empty_conformation(self):

		self.conformation.append(FFEA_script_conformation())
		self.num_conformations += 1

class FFEA_script_conformation:

	def __init__(self):
		
		self.motion_state = ""
		self.nodes = ""
		self.topology = ""
		self.surface = ""
		self.material = ""
		self.stokes = ""
		self.vdw = ""
		self.pin = ""
		self.bsites = ""

	def write_to_file(self, fout, fname, calc_kinetics):

		astr = ""
		astr += "\t\t<conformation>\n"
		astr += "\t\t\t<motion_state = %s>\n" % (self.motion_state)
		if self.motion_state != "STATIC":
			astr += "\t\t\t<topology = %s>\n" % (os.path.relpath(self.topology, os.path.dirname(os.path.abspath(fname))))
			astr += "\t\t\t<material = %s>\n" % (os.path.relpath(self.material, os.path.dirname(os.path.abspath(fname))))
			astr += "\t\t\t<stokes = %s>\n" % (os.path.relpath(self.stokes, os.path.dirname(os.path.abspath(fname))))
			astr += "\t\t\t<pin = %s>\n" % (os.path.relpath(self.pin, os.path.dirname(os.path.abspath(fname))))

		astr += "\t\t\t<nodes = %s>\n" % (os.path.relpath(self.nodes, os.path.dirname(os.path.abspath(fname))))
		astr += "\t\t\t<surface = %s>\n" % (os.path.relpath(self.surface, os.path.dirname(os.path.abspath(fname))))
		astr += "\t\t\t<vdw = %s>\n" % (os.path.relpath(self.vdw, os.path.dirname(os.path.abspath(fname))))
		if(calc_kinetics == 1 and self.bsites != ""):
			astr += "\t\t\t<binding_sites = %s>\n" % (os.path.relpath(self.bsites, os.path.dirname(os.path.abspath(fname))))
		astr += "\t\t</conformation>\n"
		fout.write(astr)

def extract_block_from_lines(title, index, lines):

	block = []
	block_index = 0
	in_block = 0
	for line in lines:
		try:
			if in_block == 0:
				if line.strip().replace("<", "").replace(">", "") == title:

					block_index += 1
					if block_index <= index:
						continue
					else:
						in_block = 1

			else:
				if line.strip().replace("<", "").replace(">", "") == "/" + title:
					break
				else:
					block.append(line)
		except:
			
			print "Error. Could not parse line " + line
			return []

	return block	
