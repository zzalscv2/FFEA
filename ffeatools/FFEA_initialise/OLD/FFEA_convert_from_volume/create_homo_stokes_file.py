import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python DAT_to_walrus_stokes_file.py [NUM NODES] [STOKES RADIUS] [OUTPUT STOKES FILE]")

num_nodes = int(sys.argv[1])
stokes_radius = float(sys.argv[2])
outfile = open(sys.argv[3], "w")
outfile.write("ffea stokes radii file\n")
s = ''
for i in xrange(num_nodes):
	s += str(stokes_radius) + "\n"

outfile.write("num_nodes " + str(num_nodes) + "\n")
outfile.write(s)
outfile.close()
