import sys, os, StringIO
from numpy import array as nparray

from FFEA_universe import *

def get_path_from_script(path, scriptdir):
	if os.path.isabs(path):
		return path
	else:
		return os.path.abspath(os.path.join(scriptdir, path))

class FFEA_script:

	def __init__(self, fname = "", fix=False):

		self.reset()

		# Start reading to test
		try:
			fin = open(fname, "r")	
		except(IOError):
			print "File " + fname  + " not found. Returning empty object..."
			self.reset()
			return

		# Get rid of all of the comments, if there are any
		fin = self.remove_all_comments(fin)
		script_lines = fin.readlines()
		fin.close()

		# Check for params block		
		if "<param>\n" not in script_lines:
			print "Error. File " + fname  + " not an FFEA script file."
			self.reset()
			return

		# Get scriptdir
		scriptdir = os.path.dirname(fname)

		# Get params
		try:
		  self.params = self.read_params_from_script_lines(script_lines, scriptdir)
		except:
			print "Error. Couldn't load <params>...</params>"
			self.reset()
			return

		if fix:
			self.params.fix_params()
		else:
			self.params.check_params()
			
		for i in range(self.params.num_blobs):
			try:
				self.blob.append(self.read_blob_from_script_lines(script_lines, scriptdir, i, self.params.num_conformations[i]))
			except:
				print "Error. Couldn't load <blob>...</blob> " + str(i)
				self.reset()
				return

		# Get springs
		try:
			self.read_springs_from_script_lines(script_lines, scriptdir)
		except:
			print "Error. Failed to load <spring>...</spring> "
			self.reset()
			return

	# # # # # # # # # # # # # # # # # # # # # #
	# we will take the comments out of iFile,
	#   write a virtual file "ffea_in" 
	#   and return its handler.
	# # # # # # # # # # # # # # # # # # # # # #
	def remove_all_comments(self, fin):
		STA = fin.readlines()
		fin.close()

		ffea_in = StringIO.StringIO()

		# and some variables to take the comments out: 
		comment = 0
		m_ini = "<!--"
		m_end = "-->"
		# Now start parsing the input file
		for txt in STA:
     
			# Strip tag wrapping 
			# line = ffea_in.readline().strip()
			line = txt.strip()
     
			# The following stuff takes care of the comments enclosed in "<!--" and "-->":
			# buf2_string = ""
			found = 0
			count = 0
			count_0 = 0
			ini = 0
			end = len(line)
			theEnd = end
     
			# remove the comments:
			while ((found != -1) and (found != len(line))):
				if (comment == 0):
					found = line.find(m_ini)
					if (found != -1):
						count += 1
						comment = 1
						ini = found
				if (comment == 1):
					found = line.find(m_end)
					if (found != -1):
						count += 2
						comment = 0
						end = found + 3
				# the line end up without closing the comment:
				if (comment == 1):
					line = line[:ini]
					break
				# we're out of the comment:
				elif (comment == 0):
					if (count == count_0 + 3):
						buf2_string = line[:ini]
						buf2_string += line[end:]
						line = buf2_string
						count_0 = count
					elif (count == count_0 + 2):
						line = line[end:]
						count_0 = count
						
			# comments removed!
			if len(line) > 0:
				ffea_in.write(line.strip() + "\n")

		ffea_in.seek(0,0)   
		return ffea_in
     
	def reset(self):
		self.params = None
		self.blob = []
		self.spring = ""

	def remove_blob(self, index=-1):
		
		self.blob.pop(index)
		self.params.num_blobs -= 1
		self.params.num_conformations.pop(index)

	def add_blob(self, blob):
		self.blob.append(blob)
		self.params.num_blobs += 1
		self.params.num_conformations.append(1)

	def read_params_from_script_lines(self, script_lines, scriptdir):

		# Get params block
		try:
			param_lines = extract_block_from_lines('param', 0, script_lines)
		except:
			param_lines = extract_block_from_lines('params', 0, script_lines)

		# Get some parameters
		params = FFEA_script_params()

		# Assign the params
		for param in param_lines:
			try:
				line = param.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].strip()
				rvalue = line.split("=")[1].strip()
				
			except:
				print "Error. Could not parse line '" + param + "'"
				continue

			
			params.assign_param(lvalue, rvalue, scriptdir = scriptdir)

		return params


	def read_blob_from_script_lines(self, script_lines, scriptdir, index, num_conformations):


		# Get relevent blob block
		blob_lines = extract_block_from_lines('blob', index, script_lines)

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
					lvalue = line.split("=")[0].strip()
					rvalue = line.split("=")[1].strip()
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
					elif lvalue == "beads":
						conformation.beads = get_path_from_script(rvalue, scriptdir)
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
				lvalue = line.split("=")[0].strip()
				rvalue = line.split("=")[1].strip()
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
		print map_lines
		blob.map_indices = []
		blob.map = []
		for line in map_lines:
			
			try:
				line = line.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].strip()
				rvalue = line.split("=")[1].strip()

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
				lvalue = line.split("=")[0].strip()
				rvalue = line.split("=")[1].strip()
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

	def read_springs_from_script_lines(self, script_lines, scriptdir):

		spring_lines = extract_block_from_lines("springs", 0, script_lines)

		if len(spring_lines) == 0:
			spring_lines = extract_block_from_lines("spring", 0, script_lines)
			if len(spring_lines) == 0:
				return

		if len(spring_lines) != 1:
			print "Error. Expected only one filename."
			return

		line = spring_lines[0]

		try:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].strip()
			rvalue = line.split("=")[1].strip()

		except(IndexError, ValueError):
			print "Error. Couldn't parse spring tag '" + line + "'"
			return

		if lvalue == "springs_fname" or lvalue == "spring_fname":
			self.spring = get_path_from_script(rvalue, scriptdir)
			
		return

	def default(self, basename):
		
		# Just in case it's not a basename...
		basename = os.path.splitext(basename)[0]

		# Default params, but change fnames
		self.params = FFEA_script_params()
		self.params.trajectory_out_fname = basename + "_trajectory.ftj"
		self.params.measurement_out_fname = basename + "_measurement.fm"
		self.params.vdw_forcefield_params = basename + ".lj"

		# Default, single blob, containing 1 conformation and no kinetics
		self.add_empty_blob()
		self.blob[-1].default(basename)

	def write_to_file(self, fname, verbose=True):

		fout = open(fname, "w")
		self.params.write_to_file(fout, fname, verbose=verbose)
		fout.write("<system>\n")
		for blob in self.blob:
			blob.write_to_file(fout, fname, self.params.calc_kinetics, self.params.calc_preComp)

		if self.spring != "":
			fout.write("\t<interactions>\n\t\t<springs>\n")
			fout.write("\t\t\t<spring_fname = %s>\n" % (os.path.relpath(self.spring, os.path.dirname(os.path.abspath(fname)))))
			fout.write("\t\t</springs>\n\t</interactions>\n")
		fout.write("</system>")
		fout.close()

	def print_details(self):
			print "traj = ", self.params.trajectory_out_fname
			print "num_blobs = ", self.params.num_blobs
			print "num_conformations = ", self.params.num_conformations
			
			for i in range(self.params.num_blobs):
				for j in range(self.params.num_conformations[i]):
					print "Node fname = ", self.blob[i].conformation[j].nodes
					

	def add_empty_blob(self):

		self.blob.append(FFEA_script_blob())
		self.params.num_blobs += 1
		self.params.num_conformations.append(1)

	# Loading other FFEA objects
	def load_node(self, bindex, cindex=0):
		return FFEA_node.FFEA_node(self.blob[bindex].conformation[cindex].nodes)

	def load_surface(self, bindex, cindex=0):
		return FFEA_surface.FFEA_surface(self.blob[bindex].conformation[cindex].surface)

	def load_topology(self, bindex, cindex=0):
		return FFEA_topology.FFEA_topology(self.blob[bindex].conformation[cindex].topology)

	def load_stokes(self, bindex, cindex=0):
		return FFEA_stokes.FFEA_stokes(self.blob[bindex].conformation[cindex].stokes)

	def load_vdw(self, bindex, cindex=0):
		return FFEA_vdw.FFEA_vdw(self.blob[bindex].conformation[cindex].vdw)

	def load_pin(self, bindex, cindex=0):
		return FFEA_pin.FFEA_pin(self.blob[bindex].conformation[cindex].pin)

	def load_material(self, bindex, cindex=0):
		return FFEA_material.FFEA_material(self.blob[bindex].conformation[cindex].material)

	def load_trajectory(self, num_frames=100000000):
		return FFEA_trajectory.FFEA_trajectory(self.params.trajectory_out_fname, num_frames_to_read = num_frames)

	def load_measurement(self, num_frames=100000000):
		return FFEA_measurement.FFEA_measurement(self.params.measurement_out_fname, num_frames_to_read = num_frames)

class FFEA_script_params():
	
	def __init__(self):

		# Must have same default values as main ffea so can load defaulted scripts
		self.restart = 0
		self.dt = 1e-14
		self.kT = 4.11e-21
		self.check = 10000
		self.num_steps = 1e11
		self.trajectory_out_fname = ""
		self.measurement_out_fname = ""
		self.vdw_forcefield_params = ""
		self.kinetics_out_fname = ""
		self.binding_site_params = ""
		self.checkpoint_in = ""
		self.checkpoint_out = ""
		self.epsilon = 0.01
		self.max_iterations_cg = 1000
		self.kappa = 1e9
		self.epsilon_0 = 1.0
		self.dielec_ext = 1
		self.calc_stokes = 1
		self.stokes_visc = 1e-3
		self.calc_vdw = 1
		self.calc_noise = 1
		self.calc_springs = 0
		self.calc_ctforces = 0
		self.inc_self_vdw = 1
		self.calc_preComp = 0
		self.calc_es = 0
		self.calc_kinetics = 0
		self.kinetics_update = 0
		self.es_update = 1
		self.es_N = nparray([-1,-1,-1])
		self.vdw_type = "steric"
		self.vdw_steric_factor = 1e-2
		self.move_into_box = 1
		self.sticky_wall_xz = 0
		self.wall_x_1 = "PBC"
		self.wall_x_2 = "PBC"
		self.wall_y_1 = "PBC"
		self.wall_y_2 = "PBC"
		self.wall_z_1 = "PBC"
		self.wall_z_2 = "PBC"
		self.es_h = 3
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
			self.check = int(float(rvalue))
		elif lvalue == "num_steps":
			self.num_steps = int(float(rvalue))
		elif lvalue == "trajectory_out_fname":
			self.trajectory_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "measurement_out_fname":
			self.measurement_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "vdw_forcefield_params":
			self.vdw_forcefield_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "kinetics_out_fname":
			self.kinetics_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "binding_site_params":
			self.binding_site_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "checkpoint_in":
			self.checkpoint_in = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "checkpoint_out":
			self.checkpoint_out = get_path_from_script(rvalue, scriptdir)
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
		elif lvalue == "calc_stokes" or lvalue == "do_stokes":
			self.calc_stokes = int(rvalue)
		elif lvalue == "calc_kinetics":
			self.calc_kinetics = int(rvalue)
		elif lvalue == "kinetics_update":
			self.kinetics_update = int(rvalue)
		elif lvalue == "stokes_visc":
			self.stokes_visc = float(rvalue)
		elif lvalue == "calc_vdw":
			self.calc_vdw = int(rvalue)
		elif lvalue == "calc_preComp":
			self.calc_preComp = int(rvalue)
		elif lvalue == "calc_springs":
			self.calc_springs = int(rvalue)
		elif lvalue == "calc_ctforces":
			self.calc_ctforces = int(rvalue)
		elif lvalue == "vdw_type":
			self.vdw_type = rvalue
		elif lvalue == "inc_self_vdw":
			self.inc_self_vdw = float(rvalue)
		elif lvalue == "vdw_steric_factor":
			self.vdw_steric_factor = float(rvalue)
		elif lvalue == "calc_noise":
			self.calc_noise = int(rvalue)
		elif lvalue == "calc_es":
			self.calc_es = int(rvalue)
		elif lvalue == "es_update":
			self.es_update = int(rvalue)
		elif lvalue == "es_N_x":
			self.es_N[0] = int(rvalue)
		elif lvalue == "es_N_y":
			self.es_N[1] = int(rvalue)
		elif lvalue == "es_N_z":
			self.es_N[2] = int(rvalue)
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
			print "Unrecognised parameter '" + lvalue + "'. Ignoring..."

	# This simply returns an error if something is incompatible. For fixes, use fix_params
	def check_params(self):
		pass

	# This function attempts to fix conflicting data (i.e. num_blobs = N but len(num_conformations) != N)
	def fix_params(self):
		if len(self.num_conformations) < self.num_blobs:
			self.num_conformations.extend([1 for i in range(self.num_blobs - len(self.num_conformations))])

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

	def write_to_file(self, fout, fname, verbose=False):
		
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
		if verbose:
			astr += "\t<dt = %5.2e>\n" % (self.dt)
			astr += "\t<kT = %5.2e>\n" % (self.kT)
			astr += "\t<check = %d>\n" % (self.check)
			astr += "\t<num_steps = %1.0e>\n" % (self.num_steps)
			astr += "\t<rng_seed = time>\n"
		astr += "\t<trajectory_out_fname = %s>\n" % (os.path.relpath(self.trajectory_out_fname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<measurement_out_fname = %s>\n" % (os.path.relpath(self.measurement_out_fname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<vdw_forcefield_params = %s>\n" % (os.path.relpath(self.vdw_forcefield_params, os.path.dirname(os.path.abspath(fname))))
		if self.kinetics_out_fname != "":
			astr += "\t<kinetics_out_fname = %s>\n" % (os.path.relpath(self.kinetics_out_fname, os.path.dirname(os.path.abspath(fname))))

		if self.binding_site_params != "":
			astr += "\t<binding_site_params = %s>\n" % (os.path.relpath(self.binding_site_params, os.path.dirname(os.path.abspath(fname))))
		if self.checkpoint_in != "":
			astr += "\t<checkpoint_in = %s>\n" % (os.path.relpath(self.checkpoint_in, os.path.dirname(os.path.abspath(fname))))
		if self.checkpoint_out != "":
			astr += "\t<checkpoint_out = %s>\n" % (os.path.relpath(self.checkpoint_out, os.path.dirname(os.path.abspath(fname))))
		if verbose:
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
			astr += "\t<calc_springs = %d>\n" % (self.calc_springs)
			astr += "\t<calc_noise = %d>\n" % (self.calc_noise)
			astr += "\t<calc_es = %d>\n" % (self.calc_es)
			astr += "\t<es_update = %d>\n" % (self.es_update)
			astr += "\t<es_N_x = %d>\n" % (self.es_N[0])
			astr += "\t<es_N_y = %d>\n" % (self.es_N[1])
			astr += "\t<es_N_z = %d>\n" % (self.es_N[2])
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
		self.scale = 1.0
		self.centroid = None
		self.rotation = None

	def write_to_file(self, fout, fname, calc_kinetics, calc_preComp):

		fout.write("\t<blob>\n")
		need_solver = 0;
		for conformation in self.conformation:
			conformation.write_to_file(fout, fname, calc_kinetics, calc_preComp)
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

	def default(self, basename):
		
		self.add_empty_conformation()
		self.conformation[-1].default(basename)

		self.solver = "CG_nomass"
		self.scale = 1e-10

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
		self.beads = ""

	def default(self, basename):

		self.motion_state = "DYNAMIC"
		self.nodes = basename + ".node"
		self.topology = basename + ".top"
		self.surface = basename + ".surf"
		self.material = basename + ".mat"
		self.stokes = basename + ".stokes"
		self.vdw = basename + ".vdw"
		self.pin = basename + ".pin"	

	def write_to_file(self, fout, fname, calc_kinetics, calc_preComp):

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
		if(calc_preComp == 1):
			astr += "\t\t\t<beads = %s>\n" % (os.path.relpath(self.beads, os.path.dirname(os.path.abspath(fname))))
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
