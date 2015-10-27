import sys, os
import FFEA_script

if len(sys.argv) != 3:
	sys.exit("Usage: python FFEA_build_simulation_environment.py [INPUT .ffea file] [Simulation root directory]")

# Get args
infname = sys.argv[1]
simdir = os.path.abspath(sys.argv[2])
scriptdir = simdir + "/scripts"
trajdir = simdir + "/traj"
measdir = simdir + "/meas"
structdir = simdir + "/structure"
resultsdir = simdir + "/results"

# Get files
if os.path.splitext(infname)[1] == ".ffea":
	ffeafname = infname
else:
	sys.exit("\nError. Unrecognised input file. Must be .ffea or .vol")

# Get script
print "Opening script..."
script = FFEA_script.FFEA_script(ffeafname)
print "done!"

# Sort files into directories
print "Making new directories and moving files..."
os.system("mkdir " + scriptdir + " " + trajdir + " " + measdir + " " + structdir + " " + resultsdir)
for b in script.blob:
	for c in b.conformation:
	
		# Nodes
		os.system("mv " + c.nodes + " " + structdir)
		c.nodes = "../structure/" + os.path.basename(os.path.abspath(c.nodes))

		# Topology
		os.system("mv " + c.topology + " " + structdir)
		c.topology = "../structure/" + os.path.basename(os.path.abspath(c.topology))

		# surface
		os.system("mv " + c.surface + " " + structdir)
		c.surface = "../structure/" + os.path.basename(os.path.abspath(c.surface))

		# material
		os.system("mv " + c.material + " " + structdir)
		c.material = "../structure/" +  os.path.basename(os.path.abspath(c.material))

		# stokes
		os.system("mv " + c.stokes + " " + structdir)
		c.stokes = "../structure/" + os.path.basename(os.path.abspath(c.stokes))

		# vdw
		os.system("mv " + c.vdw + " " + structdir)
		c.vdw = "../structure/" + os.path.basename(os.path.abspath(c.vdw))

		# pin
		os.system("mv " + c.pin + " " + structdir)
		c.pin = "../structure/" + os.path.basename(os.path.abspath(c.pin))

# Final script change
script.params.trajectory_out_fname = "../traj/" + os.path.basename(os.path.abspath(script.params.trajectory_out_fname))
script.params.measurement_out_basefname = "../meas/" + os.path.basename(os.path.abspath(script.params.measurement_out_basefname))
os.system("mv " + script.params.vdw_forcefield_params + " " + structdir)
script.params.vdw_forcefield_params = "../structure/" +  os.path.basename(os.path.abspath(script.params.vdw_forcefield_params))
print "done!"

# Finally move script
script.write_to_file(ffeafname)
os.system("mv " + ffeafname + " " + scriptdir)
print "FFEA file '" + ffeafname + "' is now in " + scriptdir
