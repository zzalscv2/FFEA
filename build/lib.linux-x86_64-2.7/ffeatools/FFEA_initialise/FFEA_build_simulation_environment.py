import sys, os
import FFEA_script

if len(sys.argv) != 3:
	sys.exit("Usage: python FFEA_build_simulation_environment.py [INPUT .ffea file] [Simulation root directory]")

# Get args
infname = sys.argv[1]
simdir = os.path.abspath(sys.argv[2])
scriptdir = simdir + "/scripts/"
trajdir = simdir + "/traj/"
measdir = simdir + "/meas/"
structdir = simdir + "/structure/"
resultsdir = simdir + "/results/"

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
		os.system("cp " + c.nodes + " " + structdir)
		c.nodes = structdir + os.path.basename(os.path.abspath(c.nodes))

		# Topology
		os.system("cp " + c.topology + " " + structdir)
		c.topology = structdir + os.path.basename(os.path.abspath(c.topology))

		# surface
		os.system("cp " + c.surface + " " + structdir)
		c.surface = structdir + os.path.basename(os.path.abspath(c.surface))

		# material
		os.system("cp " + c.material + " " + structdir)
		c.material = structdir +  os.path.basename(os.path.abspath(c.material))

		# stokes
		os.system("cp " + c.stokes + " " + structdir)
		c.stokes = structdir + os.path.basename(os.path.abspath(c.stokes))

		# vdw
		os.system("cp " + c.vdw + " " + structdir)
		c.vdw = structdir + os.path.basename(os.path.abspath(c.vdw))

		# pin
		os.system("cp " + c.pin + " " + structdir)
		c.pin = structdir + os.path.basename(os.path.abspath(c.pin))

# Final script change
script.params.trajectory_out_fname = trajdir + os.path.basename(os.path.abspath(script.params.trajectory_out_fname))
script.params.measurement_out_fname = measdir + os.path.basename(os.path.abspath(script.params.measurement_out_fname))
os.system("cp " + script.params.vdw_forcefield_params + " " + structdir)
script.params.vdw_forcefield_params = structdir +  os.path.basename(os.path.abspath(script.params.vdw_forcefield_params))
print "done!"

# Finally move script
ffeafname = scriptdir + os.path.basename(ffeafname)
script.write_to_file(ffeafname)
print "FFEA file '" + ffeafname + "' is now in " + scriptdir
