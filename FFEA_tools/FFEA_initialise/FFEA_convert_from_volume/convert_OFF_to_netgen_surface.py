import sys

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_OFF_to_netgen_surface.py [INPUT OFF FILE] [OUTPUT NETGEN SURFACE FILE]")

file_in = open(sys.argv[1],"r")
file_out = open(sys.argv[2],"w")

# read in the OFF line
line = file_in.readline()
if "OFF" not in line:
	print "Not an OFF file (first line is not OFF)"

# write the header line in the output file
file_out.write("surfacemesh\n")

# read in the line with the number of nodes, number of faces, and some mysterious number I have no idea about
line = file_in.readline()
a = line.split()
num_nodes = int(a[0])
num_faces = int(a[1])
print "num_nodes = " + str(num_nodes)
print "num_faces = " + str(num_faces)

# write the nodes section of the output file
file_out.write(str(num_nodes) + "\n")
for i in range(num_nodes):
	line = file_in.readline()
	file_out.write(line)

# write the faces section of the output file
file_out.write(str(num_faces) + "\n")
for i in range(num_faces):
	line = file_in.readline()
	a = line.split()
	n1 = str(int(a[1]) + 1)
	n2 = str(int(a[2]) + 1)
	n3 = str(int(a[3]) + 1)

	# miss out the mysterious first number
	file_out.write(n1 + " " + n2 + " " + n3 + "\n")


file_in.close()
file_out.close()
