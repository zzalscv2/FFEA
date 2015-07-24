import sys, os

def get_block_file(filein, block_title, index):
	
	block = []
	block_index = 0
	while block_index <= index:
		while filein.readline().strip().replace("<", "").replace(">", "") != block_title:
			continue

		block_index += 1

	line = filein.readline()
	while line.strip().replace("<", "").replace(">", "") != "/" +  block_title:
		block.append(line)
		line = filein.readline()
	
	return block

def get_block_lines(linesin, block_title, index):
	
	block = []
	block_index = 0
	in_block = 0
	for line in linesin:
		if in_block == 0:
			if line.strip().replace("<", "").replace(">", "") == block_title:

				block_index += 1
				if block_index <= index:
					continue
				else:
					in_block = 1

		else:
			if line.strip().replace("<", "").replace(">", "") == "/" + block_title:
				break
			else:
				block.append(line)

	return block

class FFEA_script:
	
	def __init__(self, script_fname):
		
		self.params = FFEA_script_params()
		self.params.read_params_from_script(script_fname)
		self.blob = [FFEA_script_blob(self.params.num_conformations[i]) for i in range(self.params.num_blobs)]
		for b in self.blob:
			b.read_blob_from_script(script_fname, self.blob.index(b))

	def write_to_file(self, fname):

		fout = open(fname, "w")
		self.params.write_to_file(fout)
		fout.write("<system>\n")
		for blob in self.blob:
			blob.write_to_file(fout)
		fout.write("</system>\n")
		fout.close()

class FFEA_script_params():
	
	def __init__(self):
		self.restart = 0
		self.dt = 0.0
		self.kT = 0.0
		self.check = 0
		self.num_steps = 0
		self.trajectory_out_fname = ""
		self.measurement_out_fname = ""
		self.vdw_forcefield_params = ""
		self.epsilon = 0.0
		self.max_iterations_cg = 0
		self.kappa = 0.0
		self.epsilon_0 = 0.0
		self.dielec_ext = 0.0
		self.do_stokes = 0
		self.stokes_visc = 0.0
		self.calc_vdw = 0
		self.calc_noise = 0
		self.calc_es = 0
		self.es_update = 0
		self.es_N_x = 0
		self.es_N_y = 0
		self.es_N_z = 0
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

	def write_to_file(self, fout):
		
		num_conformations_string = ""
		num_states_string = ""
		for i in range(self.num_blobs):
			num_conformations_string += str(self.num_conformations[i]) + ","
			num_states_string += str(self.num_states[i]) + ","

		num_conformations_string = num_conformations_string[0:-1]
		num_states_string = num_states_string[0:-1]

		astr = ""
		astr += "<param>\n"
		astr += "\t<restart = %d>\n" % (self.restart)
		astr += "\t<dt = %5.2e>\n" % (self.dt)
		astr += "\t<kT = %5.2e>\n" % (self.kT)
		astr += "\t<check = %d>\n" % (self.check)
		astr += "\t<num_steps = %d>\n" % (self.num_steps)
		astr += "\t<rng_seed = time>\n"
		astr += "\t<trajectory_out_fname = %s>\n" % (self.trajectory_out_fname)
		astr += "\t<measurement_out_fname = %s>\n" % (self.measurement_out_fname)
		astr += "\t<vdw_forcefield_params = %s>\n" % (self.vdw_forcefield_params)
		astr += "\t<epsilon = %5.2e>\n" % (self.epsilon)
		astr += "\t<max_iterations_cg = %d>\n" % (self.max_iterations_cg)
		astr += "\t<kappa = %5.2e>\n" % (self.kappa)
		astr += "\t<epsilon_0 = %5.2e>\n" % (self.epsilon_0)
		astr += "\t<dielec_ext = %5.2e>\n" % (self.dielec_ext)
		astr += "\t<do_stokes = %d>\n" % (self.do_stokes)
		astr += "\t<stokes_visc = %5.2e>\n" % (self.stokes_visc)
		astr += "\t<calc_vdw = %d>\n" % (self.calc_vdw)
		astr += "\t<calc_noise = %d>\n" % (self.calc_noise)
		astr += "\t<calc_es = %d>\n" % (self.calc_es)
		astr += "\t<es_update = %d>\n" % (self.es_update)
		astr += "\t<es_N_x = %d>\n" % (self.es_N_x)
		astr += "\t<es_N_y = %d>\n" % (self.es_N_y)
		astr += "\t<es_N_z = %d>\n" % (self.es_N_z)
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
		astr += "\t<num_states = (%s)>\n" % (num_states_string)
		astr += "</param>\n\n"	
		fout.write(astr)

	def read_params_from_script(self, fname):
		
		fin = open(fname, "r")
		param_lines = get_block_file(fin, "param", 0)
		for line in param_lines:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].rstrip()
			rvalue = line.split("=")[1].lstrip()

			if lvalue == "restart":
				self.restart = int(rvalue)
			elif lvalue == "dt":
				self.dt = float(rvalue)
			elif lvalue == "kT":
				self.kT = float(rvalue)
			elif lvalue == "check":
				self.check = int(rvalue)
			elif lvalue == "num_steps":
				self.num_steps = int(float(rvalue))
			elif lvalue == "trajectory_out_fname":
				self.trajectory_out_fname = rvalue
			elif lvalue == "measurement_out_fname":
				self.measurement_out_fname = rvalue
			elif lvalue == "vdw_forcefield_params":
				self.vdw_forcefield_params = rvalue
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
			elif lvalue == "do_stokes":
				self.do_stokes = int(rvalue)
			elif lvalue == "stokes_visc":
				self.stokes_visc = float(rvalue)
			elif lvalue == "calc_vdw":
				self.calc_vdw = int(rvalue)
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

		fin.close()

class FFEA_script_blob:
		
	def __init__(self, num_conformations):
		
		self.num_conformations = num_conformations
		self.conformation = [FFEA_script_conformation() for i in range(self.num_conformations)]
		self.solver = ""
		self.scale = 0.0

	def write_to_file(self, fout):

		fout.write("\t<blob>\n")
		for conformation in self.conformation:
			conformation.write_to_file(fout)

		fout.write("\t\t<solver = %s>\n" % (self.solver))
		fout.write("\t\t<scale = %5.2e>\n" % (self.scale))
		fout.write("\t</blob>\n")

	def read_blob_from_script(self, fname, blob_index):
		 
		fin = open(fname, "r")
		blob_lines = get_block_file(fin, "blob", blob_index)
		for c in self.conformation:
			c.read_conformation_from_lines(blob_lines, self.conformation.index(c))
			
		for line in blob_lines:
			if "=" not in line:
				continue

			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].rstrip()
			rvalue = line.split("=")[1].lstrip()
			if lvalue == "solver":
				self.solver = rvalue
			elif lvalue == "scale":
				self.scale = float(rvalue)

		fin.close()

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

	def write_to_file(self, fout):

		astr = ""
		astr += "\t\t<conformation>\n"
		astr += "\t\t\t<motion_state = %s>\n" % (self.motion_state)
		astr += "\t\t\t<nodes = %s>\n" % (self.nodes)
		astr += "\t\t\t<topology = %s>\n" % (self.topology)
		astr += "\t\t\t<surface = %s>\n" % (self.surface)
		astr += "\t\t\t<material = %s>\n" % (self.material)
		astr += "\t\t\t<stokes = %s>\n" % (self.stokes)
		astr += "\t\t\t<vdw = %s>\n" % (self.vdw)
		astr += "\t\t\t<pin = %s>\n" % (self.pin)
		astr += "\t\t</conformation>\n"
		fout.write(astr)

	def read_conformation_from_lines(self, blob_lines, index):

		conformation_lines = get_block_lines(blob_lines, "conformation", index)
		for line in conformation_lines:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].rstrip()
			rvalue = line.split("=")[1].lstrip()
			
			if lvalue == "motion_state":
				self.motion_state = rvalue
			elif lvalue == "nodes":
				self.nodes = rvalue
			elif lvalue == "topology":
				self.topology = rvalue
			elif lvalue == "surface":
				self.surface = rvalue
			elif lvalue == "material":
				self.material = rvalue
			elif lvalue == "stokes":
				self.stokes = rvalue
			elif lvalue == "vdw":
				self.vdw = rvalue
			elif lvalue == "pin":
				self.pin = rvalue
