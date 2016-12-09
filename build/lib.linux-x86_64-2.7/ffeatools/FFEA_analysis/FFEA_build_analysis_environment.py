import sys, os
import FFEA_script

def destroy():
	os.system("rm -r " + rootdir + "/scripts")
	os.system("rm -r " + rootdir + "/traj")
	os.system("rm -r " + rootdir + "/meas")
	os.system("rm -r " + rootdir + "/results")
	os.system("rm -r " + rootdir + "/structure")

if len(sys.argv) != 3:
	sys.exit("Usage: python FFEA_build_analysis_environment.py [INPUT .ffea script] [Analysis root directory]")

# Get args and expand
infname = os.path.abspath(sys.argv[1])
oldscriptdir = os.path.dirname(infname)
rootdir = os.path.abspath(sys.argv[2])
trajdir = rootdir + "/traj/"
measdir = rootdir + "/meas/"
structdir = rootdir + "/structure/"
scriptdir = rootdir + "/scripts/"
resultdir = rootdir + "/results/"

# Get a script object
script = FFEA_script.FFEA_script(infname)

# Build directory structure
os.system("mkdir " + trajdir + " " + measdir + " " + scriptdir + " " + resultdir + " " + structdir)

# Now, for all relevent things in the script file, copy them to the new directories and rename them in the script

# Traj
os.system("cp " + script.params.trajectory_out_fname + " " + trajdir)
script.params.trajectory_out_fname = "../traj/" + os.path.basename(script.params.trajectory_out_fname)

# Meas
for i in range(len(script.params.measurement_out_fname)):
	os.system("cp " + script.params.measurement_out_fname[i] + " " + measdir)
	script.params.measurement_out_fname[i] = "../meas/" + os.path.basename(script.params.measurement_out_fname[i])

script.params.measurement_out_basefname = os.path.splitext(script.params.measurement_out_fname[-1])[0][0:-6] + os.path.splitext(script.params.measurement_out_fname[-1])[1]

# LJ
os.system("cp " + script.params.vdw_forcefield_params + " " + structdir)
script.params.vdw_forcefield_params = "../structure/" +  os.path.basename(script.params.vdw_forcefield_params)

# Structure
for b in script.blob:
	for c in b.conformation:

		# Nodes
		os.system("cp " + c.nodes + " " + structdir)
		c.nodes = "../structure/" + os.path.basename(c.nodes)

		# Topology
		os.system("cp " + c.topology + " " + structdir)
		c.topology = "../structure/" + os.path.basename(c.topology)

		# surface
		os.system("cp " + c.surface + " " + structdir)
		c.surface = "../structure/" + os.path.basename(c.surface)

		# material
		os.system("cp " + c.material + " " + structdir)
		c.material = "../structure/" +  os.path.basename(c.material)

		# stokes
		os.system("cp " + c.stokes + " " + structdir)
		c.stokes = "../structure/" + os.path.basename(c.stokes)

		# vdw
		os.system("cp " + c.vdw + " " + structdir)
		c.vdw = "../structure/" + os.path.basename(c.vdw)

		# pin
		os.system("cp " + c.pin + " " + structdir)
		c.pin = "../structure/" + os.path.basename(c.pin)

# Script
outfname = scriptdir + os.path.basename(infname)
script.write_to_file(outfname)
