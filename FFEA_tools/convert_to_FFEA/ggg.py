import os, sys

if len(sys.argv) != 6 and len(sys.argv) != 2:
	sys.exit("Usage: python convert_all_to_walrus.py [INPUT FILE] (OPTIONAL){[OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE] [OUTPUT SURFACE FILE] [SCALING FACTOR]}")

path = os.path.dirname(sys.argv[0])
if path == '':
	path = '.'

print "path of running script: " + path

inputfile = sys.argv[1]
print "Input file = " + inputfile
root, ext = os.path.splitext(inputfile)
base = os.path.basename(root)
print "root = " + root
print "extension = " + ext
print "base = " + base

if len(sys.argv) == 6:
	node = sys.argv[2]
	top = sys.argv[3]
	surf = sys.argv[4]
	scaling_factor = sys.argv[5]

if len(sys.argv) == 2:
	print "Only input file name given. Generating output file names from the extensionless basename."
	node = base + ".node"	
	top = base + ".top"	
	surf = base + ".surf"	
	scaling_factor = "1"

print "Script will run with the following parameters:"
print "node file = " + node
print "topology file = " + top
print "surface file = " + surf
print "scaling factor = " + scaling_factor

#print "Determining file type from extension:"
#if ext == ".vol":
#	print "VOL extension = Netgen mesh output file. Converting to walrus format..."
#	os.system("python " + path + "/netgen_to_walrus.py " + inputfile + " " + node + " " + top + "\n")
#elif ext == ".off":
#	print "OFF extension = DEC Object File Format."
#	print "Warning: This is a surface mesh, and needs to be meshed."
#	print "Meshing with '~/tetgen1.4.3/tetgen -a -Q" + inputfile + "'..."
#	os.system("~/tetgen1.4.3/tetgen -a -Q " + inputfile)
#	print "Now converting tetgen output to walrus format..."
#	os.system("python " + path + "/tetgen_to_walrus.py " + root + ".1.node " + root + ".1.ele " + node + " " + top + "\n")
#else:
#	print "Error: File extension not recognised. Only '.vol' and '.off' are supported."
#	sys.exit(0)

# convert 4 node tets to 10 point tets
os.system("python " + path + "/convert_tetrahedra_linear_to_quadratic.py " + node + " " + top + "\n")

# extract the surface
os.system("python " + path + "/extract_surface_from_topology.py " + top + " " + surf + "\n")

# move all surface nodes to the top of the list
os.system("python " + path + "/put_surface_nodes_first.py " + node + " " + top + " " + surf + "\n")

# move all surface elements to the top of the list
os.system("python " + path + "/put_surface_elements_first.py " + node + " " + top + " " + surf + "\n")

# cycle surface indices so that all normals point out of the mesh
os.system("python " + path + "/make_all_normals_point_outwards.py " + node + " " + top + " " + surf + "\n")

# calculate the element connectivity (for use with van der waals calculations)
#os.system("python " + path + "/calculate_element_connectivity.py " + top + " " + vdw + "\n")

# scale the cube by the desired factor
os.system("python " + path + "/scale.py " + node + " " + node + " " + scaling_factor + "\n")

print "All done"
