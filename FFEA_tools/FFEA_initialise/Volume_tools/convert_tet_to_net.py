import sys, os
import FFEA_node, FFEA_topology, FFEA_surface

def error(m = ""):
	print "\n\tERROR: " + m + "\n"
	sys.exit()

def usage():
	print "Usage: python " + os.path.basename(os.path.abspath(sys.argv[0]) + "-i <input_file> -o <output_file>")

if len(sys.argv) == 1 or (len(sys.argv) - 1) % 2 != 0:
	usage()
	sys.exit()

# Get args
input_fnames = []
output_fname = ""
for i in range(1, len(sys.argv), 2):
	if sys.argv[i] == "-i":
		input_fnames.append(sys.argv[i + 1])
	elif sys.argv[i] == "-o":
		output_fname = sys.argv[i + 1]

# Read files
top = None
surf = None
node = None

for f in input_fnames:
	base, ext = os.path.splitext(f)

	if ext == ".ele":
		top = FFEA_topology.FFEA_topology(f)

		# Additionally
		if output_fname == "":
			output_fname = base + ".vol"

	elif ext == ".node":
		node = FFEA_node.FFEA_node(f)
	elif ext == ".face":
		surf = FFEA_surface.FFEA_surface(f)

if top == None or surf == None or node == None:
	error("A topology (.ele), surface (.face) and node (.node) are all required to create a NETGEN .vol file.")

# Now, write them all to file
fout = open(output_fname, "w")
fout.write("mesh3d\ndimension\n3\ngeomtype\n11\n\n")
fout.close()

surf.write_to_file(output_fname)
top.write_to_file(output_fname)
node.write_to_file(output_fname)

fout = open(output_fname, "a")
fout.write("endmesh\n")
fout.close()
