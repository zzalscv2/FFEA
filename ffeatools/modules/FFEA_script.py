# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

import sys, os
if sys.version_info[0] < 3:
	from StringIO import StringIO
else:
	from io import StringIO

from numpy import array as nparray

from FFEA_universe import *

def get_path_from_script(path, scriptdir):
	if os.path.isabs(path):
		return path
	else:
		return os.path.abspath(os.path.join(scriptdir, path))

class FFEA_script:

	def __init__(self, fname = "", fix=True):

		self.reset()
		sys.stdout.write("Loading FFEA script file...")

		# Return empty object if fname not initialised
		if fname == "":
			self.valid = True
			sys.stdout.write("Empty script object initialised.\n")
			return

		# Start reading to test
		try:
			fin = open(fname, "r")	
		except(IOError):
			raise

		# Get rid of all of the comments, if there are any
		fin = self.remove_all_comments(fin)
		script_lines = fin.readlines()
		fin.close()

		# Check for params block		
		if "<param>\n" not in script_lines:
			print("Error. File " + fname  + " not an FFEA script file.")
			self.reset()
			return

		# Get scriptdir
		scriptdir = os.path.dirname(fname)

		# Get params
		try:
		  self.params = self.read_params_from_script_lines(script_lines, scriptdir)
		except:
			print("Error. Couldn't load <params>...</params>")
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
				print("Error. Couldn't load <blob>...</blob> " + str(i))
				self.reset()
				return
        
        # load rods
		try:
			for i in range(self.params.num_rods):
				self.rod.append(self.read_rod_from_script_lines(script_lines, scriptdir, i, ))
		#except AttributeError:
		#	pass # no rods exist
		except:
				print("Error. Couldn't load <rod>...</rod> " + str(i))
				self.reset()
				raise # why do all these functions return instead of raising...?
        
		# Get springs
		try:
			self.read_springs_from_script_lines(script_lines, scriptdir)
		except:
			print("Error. Failed to load <spring>...</spring> ")
			self.reset()
			return

		# And forces
		try:
			self.read_ctforces_from_script_lines(script_lines, scriptdir)
		except:
			print("Error. Failed to load <ctforces>...</ctforces> ")
			self.reset()
			return

		# Get precomp-stuff
		self.precomp = FFEA_script_precomp()
		try:
			self.precomp.read_precomp_from_script_lines(script_lines, scriptdir)
		except:
			print("Error. Failed to load <precomp>...</precomp> ")
			self.reset()
			return

		if self.precomp.types == [] or self.precomp.types[0] == "":
			self.precomp = None

		self.valid = True
		self.empty = False
		sys.stdout.write("done!\n")

	# # # # # # # # # # # # # # # # # # # # # #
	# we will take the comments out of iFile,
	#   write a virtual file "ffea_in" 
	#   and return its handler.
	# # # # # # # # # # # # # # # # # # # # # #
	def remove_all_comments(self, fin):
		STA = fin.readlines()
		fin.close()

		ffea_in = StringIO()

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
				print("Error. Could not parse line '" + param + "'")
				continue

			
			params.assign_param(lvalue, rvalue, scriptdir = scriptdir)

		return params


	def read_blob_from_script_lines(self, script_lines, scriptdir, index, num_conformations):

		# Get relevent blob block
		blob_lines = extract_block_from_lines('blob', index, script_lines)

		# Get some parameters
		blob = FFEA_script_blob()

		# Get conformations first
		enforce_conf_blocks = False
		if (num_conformations > 1): enforce_conf_blocks = True
		for i in range(num_conformations):

			# Get a conformation
			conformation = FFEA_script_conformation()
			
			conformation_lines = extract_block_from_lines('conformation', i, blob_lines)
			read_blob_as_conf = False
			if len(conformation_lines) == 0:
				if enforce_conf_blocks == True: 
					print("Error. Couldn't parse conformation tag '" + line + "'")
					return None
				else:
					conformation_lines = blob_lines[:]
					read_blob_as_conf = True

			for line in conformation_lines:
				try:
					line = line.strip().replace("<", "").replace(">", "")
					lvalue = line.split("=")[0].strip()
					rvalue = line.split("=")[1].strip()
				except:
					print("Error. Could not parse conformation tag '" + line + "'")
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
					elif lvalue == "skeleton":
						conformation.skeleton = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "binding_sites":
						conformation.bsites = get_path_from_script(rvalue, scriptdir)
					elif lvalue == "beads":
						conformation.beads = get_path_from_script(rvalue, scriptdir)
					else:
						if ((read_blob_as_conf == True) and ( (lvalue == "centroid") or \
							(lvalue == "rotation") or (lvalue == "solver") or (lvalue == "scale") or \
							(lvalue == "translate") or (lvalue == "velocity") or (lvalue == "velocity") or \
                   	(lvalue == "calc_compress") or (lvalue == "compres"))): continue
						else: 
							print("Unrecognised conformation tag '" + line + "'. Ignoring...")
							continue

				except(IndexError, ValueError):
					print("Error. Couldn't parse conformation tag '" + line + "'")
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
				print("Error. Could not parse blob tag '" + line + "'")
				return None

			try:
				if lvalue == "states":
					blob.states = rvalue
				elif lvalue == "rates":
					blob.rates = rvalue
				else:
					continue

			except(IndexError, ValueError):
				print("Error. Couldn't parse blob tag '" + line + "'")
				return None

		# Now maps
		map_lines = extract_block_from_lines('maps', 0, kinetic_lines)
		blob.map_indices = []
		blob.map = []
		for line in map_lines:
			
			try:
				line = line.strip().replace("<", "").replace(">", "")
				lvalue = line.split("=")[0].strip()
				rvalue = line.split("=")[1].strip()

			except:
				print("Error. Could not parse blob tag '" + line + "'")
				return None

			try:
				blob.map_indices.append([int(lvalue.strip().split("(")[1].replace(")","").split(",")[i]) for i in range(2)])
				blob.map_indices
				blob.map.append(rvalue)

			except:
				print("Error. Could not parse blob tag '" + line + "'")
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
				print("Error. Could not parse blob tag '" + line + "'")
				return None

			try:
				if lvalue == "solver":
					if(rvalue.strip() != ""):
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
				print("Error. Couldn't parse blob tag '" + line + "'")
				return None

		return blob
    
	def read_rod_from_script_lines(self, script_lines, scriptdir, index):
         # Get the starting and ending line numbers for the <rod> block
         rod_index = 0
         for line_no in range(len(script_lines)):
            if "<rod>" in script_lines[line_no]:
                if rod_index == index:
                    line_start = line_no
            if "</rod>" in script_lines[line_no]:
                if rod_index == index:
                    line_end = line_no
                rod_index +=1
         
         # Can't avoid it now, just gotta do the FFEA_script thing
         # (aka: write a bunch of spaghetti code instead of just using an XML parser)
         # First, parse thing
         rod_lines = script_lines[line_start+1:line_end]
         for line_no in range(len(rod_lines)):
             rod_lines[line_no] = rod_lines[line_no].replace("\t", "")
             rod_lines[line_no] = rod_lines[line_no].replace("\n", "")
             rod_lines[line_no] = rod_lines[line_no].replace("<", "")
             rod_lines[line_no] = rod_lines[line_no].replace(">", "")
        
         #absolute vs relative paths (don't begin a relative path with a /)
         if scriptdir == "":
             separator = ""
         else:
             separator = "/"
        
         # Read in rod
         for line in rod_lines:
             if line.split("=")[0].strip() == "input":
                 in_path = scriptdir+separator+line.split("=")[1].strip()
            
             if line.split("=")[0].strip() == "output":
                 out_path = scriptdir+separator+line.split("=")[1].strip()
                 
         try:
             rod = FFEA_rod.FFEA_rod(out_path)
         except IOError:
             print("Couldn't find "+out_path)
             rod = FFEA_rod.FFEA_rod(in_path)

          #scaling, translation, rotation etc go here
         return rod
             
         
         

	def read_springs_from_script_lines(self, script_lines, scriptdir):

		spring_lines = extract_block_from_lines("springs", 0, script_lines)

		if len(spring_lines) == 0:
			spring_lines = extract_block_from_lines("spring", 0, script_lines)
			if len(spring_lines) == 0:
				return

		if len(spring_lines) != 1:
			print("Error. Expected only one filename.")
			return

		line = spring_lines[0]

		try:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].strip()
			rvalue = line.split("=")[1].strip()

		except(IndexError, ValueError):
			print("Error. Couldn't parse spring tag '" + line + "'")
			return

		if lvalue == "springs_fname" or lvalue == "spring_fname":
			self.spring = get_path_from_script(rvalue, scriptdir)
			print(self.spring)
			
		return

	def read_ctforces_from_script_lines(self, script_lines, scriptdir):

		ctforce_lines = extract_block_from_lines("ctforces", 0, script_lines)

		if len(ctforce_lines) == 0:
			ctforce_lines = extract_block_from_lines("ctforce", 0, script_lines)
			if len(ctforce_lines) == 0:
				return

		if len(ctforce_lines) != 1:
			print("Error. Expected only one filename.")
			return

		line = ctforce_lines[0]

		try:
			line = line.strip().replace("<", "").replace(">", "")
			lvalue = line.split("=")[0].strip()
			rvalue = line.split("=")[1].strip()

		except(IndexError, ValueError):
			print("Error. Couldn't parse ctforce tag '" + line + "'")
			return

		if lvalue == "ctforce_fname" or lvalue == "ctforces_fname":
			self.ctforces = get_path_from_script(rvalue, scriptdir)
			
		return



	def default(self, basename, checkBasename=True):
		
		if checkBasename:
			# Just in case it's not a basename...
			basename = os.path.splitext(basename)[0]

		# Default params, but change fnames
		self.params = FFEA_script_params()
		self.params.trajectory_out_fname = basename + ".ftj"
		self.params.measurement_out_fname = basename + ".fm"
		self.params.vdw_forcefield_params = basename + ".lj"

		# Default, single blob, containing 1 conformation and no kinetics
		self.add_empty_blob()
		self.blob[-1].default(basename)

	def write_to_file(self, fname, verbose=False):

		fout = open(fname, "w")
		self.params.write_to_file(fout, fname, verbose=verbose)
		fout.write("<system>\n")
		for blob in self.blob:
			blob.write_to_file(fout, fname, self.params.calc_kinetics, self.params.calc_preComp, verbose = verbose)


		# I think we should put all of this in an 'interactions' object. Then we can have 'if interactions == None' or something. I can't though, no time, writing my thesis :D
		if self.spring != "" or self.precomp != None or self.ctforces != "":
			fout.write("\t<interactions>\n")
			if self.spring != "":
				fout.write("\t\t<springs>\n")
				fout.write("\t\t\t<spring_fname = %s>\n" % (os.path.relpath(self.spring, os.path.dirname(os.path.abspath(fname)))))
				fout.write("\t\t</springs>\n")
			if self.ctforces != "":
				fout.write("\t\t<ctforces>\n")
				fout.write("\t\t\t<ctforces_fname = %s>\n" % (os.path.relpath(self.ctforces, os.path.dirname(os.path.abspath(fname)))))
				fout.write("\t\t</ctforces>\n")

			if self.precomp != None:
				self.precomp.write_to_file(fout)

			fout.write("\t</interactions>\n")
		fout.write("</system>\n")
		fout.close()

	def print_details(self):
			print("traj = ", self.params.trajectory_out_fname)
			print("num_blobs = ", self.params.num_blobs)
			print("num_conformations = ", self.params.num_conformations)
			
			for i in range(self.params.num_blobs):
				for j in range(self.params.num_conformations[i]):
					print("Node fname = ", self.blob[i].conformation[j].nodes)
					

	def add_empty_blob(self):

		self.blob.append(FFEA_script_blob())
		self.params.num_blobs += 1
		self.params.num_conformations.append(1)

	def reset(self):
		self.valid = False
		self.empty = True
		self.params = None
		self.blob = []
		self.rod = []
		self.spring = ""
		self.ctforces = ""
		self.precomp = None

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

	def load_skeleton(self, bindex, cindex=0):
		return FFEA_skeleton.FFEA_skeleton(self.blob[bindex].conformation[cindex].skeleton)

	def load_material(self, bindex, cindex=0):
		return FFEA_material.FFEA_material(self.blob[bindex].conformation[cindex].material)

	def load_trajectory(self, num_frames=100000000, start=0, frame_rate = 1):
		return FFEA_trajectory.FFEA_trajectory(self.params.trajectory_out_fname, num_frames_to_read = num_frames, start=start, frame_rate = frame_rate)

	def load_measurement(self, num_frames=100000000):
		return FFEA_measurement.FFEA_measurement(self.params.measurement_out_fname, num_frames_to_read = num_frames)

	def load_lj(self):
		return FFEA_lj.FFEA_lj(self.params.vdw_forcefield_params)

	def load_ctforces(self):
		return FFEA_ctforces.FFEA_ctforces(self.ctforces)

class FFEA_script_params():
	
	def __init__(self):

		# Must have same default values as main ffea so can load defaulted scripts
		self.restart = 0
		self.dt = 1e-14
		self.kT = 4.11e-21
		self.rng_seed = "" # string
		self.kt = self.kT
		self.check = 10000
		self.num_steps = 1e11
		self.trajectory_out_fname = ""
		self.measurement_out_fname = ""
		self.vdw_forcefield_params = ""
		self.kinetics_out_fname = ""
		self.binding_site_params = ""
		self.beads_out_fname = ""
		self.checkpoint_in = ""
		self.checkpoint_out = ""
		self.epsilon = 0.01
		self.max_iterations_cg = 1000
		self.kappa = 1e9
		self.epsilon_0 = 1.0
		self.dielec_ext = 1
		self.calc_stokes = 1
		self.stokes_visc = 1e-3
		self.calc_ssint = 1
		self.calc_noise = 1
		self.calc_springs = 0
		self.calc_ctforces = 0
		self.inc_self_vdw = 1
		self.ssint_cutoff = nparray([3.0e-9, 3.0e-9, 3.0e-9])
		self.calc_preComp = 0
		self.calc_es = 0
		self.calc_kinetics = 0
		self.kinetics_update = 0
		self.es_update = 1
		self.es_N = nparray([-1,-1,-1])
		self.ssint_type = "ljsteric"
		self.steric_factor = 1
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
			self.rng_seed = rvalue
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
		elif lvalue == "vdw_forcefield_params" or lvalue == "vdw_in_fname" or lvalue == "lj_params" or lvalue == "ljparams":
			self.vdw_forcefield_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "kinetics_out_fname":
			self.kinetics_out_fname = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "binding_site_params":
			self.binding_site_params = get_path_from_script(rvalue, scriptdir)
		elif lvalue == "beads_out_fname":
			self.beads_out_fname = get_path_from_script(rvalue, scriptdir)
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
		elif lvalue == "calc_ssint":
			self.calc_ssint = int(rvalue)
		elif lvalue == "calc_preComp":
			self.calc_preComp = int(rvalue)
		elif lvalue == "calc_springs":
			self.calc_springs = int(rvalue)
		elif lvalue == "calc_ctforces":
			self.calc_ctforces = int(rvalue)
		elif lvalue == "ssint_type":
			self.ssint_type = rvalue
		elif lvalue == "inc_self_vdw":
			self.inc_self_vdw = int(rvalue)
		elif lvalue == "ssint_cutoff":
			self.ssint_cutoff = nparray([float(rvalue) for i in range(3)])
		elif lvalue == "ssint_cutoff_x":
			self.ssint_cutoff[0] = float(rvalue)
		elif lvalue == "ssint_cutoff_y":
			self.ssint_cutoff[1] = float(rvalue)
		elif lvalue == "ssint_cutoff_z":
			self.ssint_cutoff[2] = float(rvalue)
		elif lvalue == "steric_factor":
			self.steric_factor = float(rvalue)
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
		elif lvalue == "num_rods":
			self.num_rods = int(rvalue)
		elif lvalue == "num_conformations":
			self.num_conformations = [int(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
		elif lvalue == "num_states":
			self.num_states = [int(r) for r in rvalue.replace("(", "").replace(")", "").split(",")]
		else:
			print("Unrecognised parameter '" + lvalue + "'. Ignoring...")

	# This simply returns an error if something is incompatible. For fixes, use fix_params
	def check_params(self):
		pass

	# This function attempts to fix conflicting data (i.e. num_blobs = N but len(num_conformations) != N)
	def fix_params(self):
		if len(self.num_conformations) < self.num_blobs:
			self.num_conformations.extend([1 for i in range(self.num_blobs - len(self.num_conformations))])
         
		try:
			self.num_rods
		except AttributeError:
			self.num_rods = 0

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
		
		write_num_conformations_string = False 
		num_conformations_string = ""
		num_states_string = ""
		for i in range(self.num_blobs):
			num_conformations_string += str(self.num_conformations[i]) + ","
			if (self.num_conformations[i]) > 1: write_num_conformations_string = True
			#num_states_string += str(self.num_states[i]) + ","

		num_conformations_string = num_conformations_string[0:-1]
		#num_states_string = num_states_string[0:-1]
		
		astr = ""
		astr += "<param>\n"
		astr += "\t<restart = %d>\n" % (self.restart)
		if verbose:
			astr += "\t<dt = %6.2e>\n" % (self.dt)
			astr += "\t<kT = %6.2e>\n" % (self.kT)
			astr += "\t<check = %d>\n" % (self.check)
			astr += "\t<num_steps = %6.2e>\n" % (self.num_steps)
			astr += "\t<rng_seed = %s>\n" % (self.rng_seed)
		astr += "\t<trajectory_out_fname = %s>\n" % (os.path.relpath(self.trajectory_out_fname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<measurement_out_fname = %s>\n" % (os.path.relpath(self.measurement_out_fname, os.path.dirname(os.path.abspath(fname))))
		astr += "\t<lj_params = %s>\n" % (os.path.relpath(self.vdw_forcefield_params, os.path.dirname(os.path.abspath(fname))))
		if self.kinetics_out_fname != "":
			astr += "\t<kinetics_out_fname = %s>\n" % (os.path.relpath(self.kinetics_out_fname, os.path.dirname(os.path.abspath(fname))))

		if self.binding_site_params != "":
			astr += "\t<binding_site_params = %s>\n" % (os.path.relpath(self.binding_site_params, os.path.dirname(os.path.abspath(fname))))
		if self.checkpoint_in != "":
			astr += "\t<checkpoint_in = %s>\n" % (os.path.relpath(self.checkpoint_in, os.path.dirname(os.path.abspath(fname))))
		if self.checkpoint_out != "":
			astr += "\t<checkpoint_out = %s>\n" % (os.path.relpath(self.checkpoint_out, os.path.dirname(os.path.abspath(fname))))
		if self.beads_out_fname != "":
			astr += "\t<beads_out_fname = %s>\n" % (os.path.relpath(self.beads_out_fname, os.path.dirname(os.path.abspath(fname))))

		if (self.calc_es == 1):
			astr += "\t<calc_es = %d>\n" % (self.calc_es)
			astr += "\t<epsilon_0 = %6.2e>\n" % (self.epsilon_0)
			astr += "\t<dielec_ext = %6.2e>\n" % (self.dielec_ext)
			astr += "\t<kappa = %6.2e>\n" % (self.kappa)
			astr += "\t<es_h = %d>\n" % (self.es_h)
		if verbose:
			astr += "\t<max_iterations_cg = %d>\n" % (self.max_iterations_cg)
			astr += "\t<epsilon = %6.2e>\n" % (self.epsilon)
		if (self.calc_stokes == 1):
			astr += "\t<calc_stokes = %d>\n" % (self.calc_stokes)
			astr += "\t<stokes_visc = %6.2e>\n" % (self.stokes_visc)
		if (self.calc_ssint == 1):
			astr += "\t<calc_ssint = %d>\n" % (self.calc_ssint)
			astr += "\t<ssint_type = %s>\n" % (self.ssint_type)
			if self.ssint_type == "steric" or self.ssint_type == "ljsteric":
				astr += "\t<steric_factor = %6.2e>\n" % (self.steric_factor)
			astr += "\t<ssint_cutoff_x = %6.2e>\n" % (self.ssint_cutoff[0])
			astr += "\t<ssint_cutoff_y = %6.2e>\n" % (self.ssint_cutoff[1])
			astr += "\t<ssint_cutoff_z = %6.2e>\n" % (self.ssint_cutoff[2])
			astr += "\t<inc_self_vdw = %d>\n" % (self.inc_self_vdw)
		if (self.calc_springs == 1):
			astr += "\t<calc_springs = %d>\n" % (self.calc_springs)
		if (self.calc_noise == 0):
			astr += "\t<calc_noise = %d>\n" % (self.calc_noise)
		if (self.calc_preComp == 1):
			astr += "\t<calc_preComp = %d>\n" % (self.calc_preComp)
		if ((self.es_N[0] != -1) or (self.es_N[1] != -1) or (self.es_N[2] != -1)):
			astr += "\t<es_N_x = %d>\n" % (self.es_N[0])
			astr += "\t<es_N_y = %d>\n" % (self.es_N[1])
			astr += "\t<es_N_z = %d>\n" % (self.es_N[2])
			astr += "\t<move_into_box = %d>\n" % (self.move_into_box)
			# astr += "\t<sticky_wall_xz = %d>\n" % (self.sticky_wall_xz)
			# astr += "\t<wall_x_1 = %s>\n" % (self.wall_x_1)
			# astr += "\t<wall_x_2 = %s>\n" % (self.wall_x_2)
			# astr += "\t<wall_y_1 = %s>\n" % (self.wall_y_1)
			# astr += "\t<wall_y_2 = %s>\n" % (self.wall_y_2)
			# astr += "\t<wall_z_1 = %s>\n" % (self.wall_z_1)
			# astr += "\t<wall_z_2 = %s>\n" % (self.wall_z_2)
		if verbose:
			astr += "\t<es_update = %d>\n" % (self.es_update)
			
		astr += "\t<num_blobs = %d>\n" % (self.num_blobs)
		if (write_num_conformations_string):
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
		self.solver = "CG_nomass"
		self.scale = 1.0
		self.centroid = []
		self.rotation = []

	def write_to_file(self, fout, fname, calc_kinetics, calc_preComp, verbose = False):

		fout.write("\t<blob>\n")
		need_solver = 0;
		need_conformations = False
		if (self.num_conformations > 1): need_conformations = True
		for conformation in self.conformation:
			conformation.write_to_file(fout, fname, calc_kinetics, calc_preComp, need_conformations, verbose = verbose)
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

		if need_solver == 1 and verbose:
			fout.write("\t\t<solver = %s>\n" % (self.solver))

		fout.write("\t\t<scale = %6.2e>\n" % (self.scale))
		if len(self.rotation) != 0:
			fout.write("\t\t<rotation = (")
			for i in range(len(self.rotation)):
				fout.write("%6.2f" % (self.rotation[i]))
				if i != len(self.rotation) - 1:
					fout.write(",")

			fout.write(")>\n")
		if len(self.centroid) != 0:
			fout.write("\t\t<centroid = (")
			for i in range(len(self.centroid)):
				fout.write("%6.2f" % (self.centroid[i]))
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
		self.skeleton = ""

	def default(self, basename):

		self.motion_state = "DYNAMIC"
		self.nodes = basename + ".node"
		self.topology = basename + ".top"
		self.surface = basename + ".surf"
		self.material = basename + ".mat"
		self.stokes = basename + ".stokes"
		self.vdw = basename + ".vdw"
		self.pin = basename + ".pin"
		self.skeleton = basename + ".skel"

	def write_to_file(self, fout, fname, calc_kinetics, calc_preComp, need_conformations=True, verbose = False):

		astr = ""
		if (need_conformations):
			astr += "\t\t<conformation>\n"
			tabs = "\t\t\t"
		else:
			tabs = "\t\t"

		astr += tabs + "<motion_state = %s>\n" % (self.motion_state)
		if self.motion_state != "STATIC":
			astr += tabs + "<topology = %s>\n" % (os.path.relpath(self.topology, os.path.dirname(os.path.abspath(fname))))
			astr += tabs + "<material = %s>\n" % (os.path.relpath(self.material, os.path.dirname(os.path.abspath(fname))))
			astr += tabs + "<stokes = %s>\n" % (os.path.relpath(self.stokes, os.path.dirname(os.path.abspath(fname))))
			if verbose and self.pin != "":
				astr += tabs + "<pin = %s>\n" % (os.path.relpath(self.pin, os.path.dirname(os.path.abspath(fname))))

		astr += tabs + "<nodes = %s>\n" % (os.path.relpath(self.nodes, os.path.dirname(os.path.abspath(fname))))
		astr += tabs + "<surface = %s>\n" % (os.path.relpath(self.surface, os.path.dirname(os.path.abspath(fname))))
		astr += tabs + "<vdw = %s>\n" % (os.path.relpath(self.vdw, os.path.dirname(os.path.abspath(fname))))

		if verbose and self.skeleton != "":
			astr += tabs + "<skeleton = %s>\n" % (os.path.relpath(self.skeleton, os.path.dirname(os.path.abspath(fname))))
	
		if (self.beads != ""):
			astr += tabs + "<beads = %s>\n" % (os.path.relpath(self.beads, os.path.dirname(os.path.abspath(fname))))


		if(calc_kinetics == 1 and self.bsites != ""):
			astr += tabs + "<binding_sites = %s>\n" % (os.path.relpath(self.bsites, os.path.dirname(os.path.abspath(fname))))

		if (need_conformations):
			astr += "\t\t</conformation>\n"

		fout.write(astr)


class FFEA_script_precomp():
	
	def __init__(self):

		# Must have same default values as main ffea so can load defaulted scripts
		self.types = []
		self.inputData = 1
		self.folder = ""
		self.dist_to_m = 1
		self.E_to_J = 1


	def read_precomp_from_script_lines(self, script_lines, scriptdir):

		precomp_lines = extract_block_from_lines("precomp", 0, script_lines)

		if len(precomp_lines) == 0:
			return

		done = 0
		for l in precomp_lines:
			try:
				l = l.strip().replace("<", "").replace(">", "")
				lvalue = l.split("=")[0].strip()
				rvalue = l.split("=")[1].strip()
	
			except(IndexError, ValueError):
				print("Error. Couldn't parse precomp tag '" + l + "'")
				return

			if lvalue == "types":
				self.types = [r.strip() for r in rvalue.replace("(", "").replace(")", "").split(",")]
				done += 1
			if lvalue == "inputData":
				self.inputData = int(rvalue)
				done += 1
			if lvalue == "folder":
				self.folder = rvalue
				done += 1
			if lvalue == "dist_to_m":
				self.dist_to_m = float(rvalue)
				done += 1
			if lvalue == "E_to_J":
				self.E_to_J = float(rvalue)
				done += 1
			
		if done != 5 : 
			print("Error. Could not parse all the PreComp fields.")
			return

		return

	def write_to_file(self, fout):
		fout.write("\t\t<precomp>\n")

		txt_types = "("
		for t in self.types:
			txt_types += t + ", "
		txt_types = txt_types[:-2] + ")"
		fout.write("\t\t\t<types = %s>\n" % (txt_types))

		fout.write("\t\t\t<inputData = %d>\n" % (self.inputData))
		fout.write("\t\t\t<folder = %s>\n" % (self.folder))
		fout.write("\t\t\t<dist_to_m = %e>\n" % (self.dist_to_m))
		fout.write("\t\t\t<E_to_J = %e>\n" % (self.E_to_J))

		fout.write("\t\t</precomp>\n")
		return


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
			
			print("Error. Could not parse line " + line)
			return []

	return block	

def check_if_block_in_lines(title, index, lines):

	block_index = 0
	block_in_lines = False

	for line in lines:
		try:
			if in_block == 0:
				if line.strip().replace("<", "").replace(">", "") == title:

					block_index += 1
					if block_index <= index:
						continue
					else:
						block_in_lines = True
						break

		except:
			
			print("Error. Could not parse line " + line)
			break

	return block_in_lines
