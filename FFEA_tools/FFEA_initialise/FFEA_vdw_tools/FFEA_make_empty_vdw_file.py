import sys, os
import FFEA_vdw

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT num_faces] [OUTPUT .vdw fname]")

# Get args
num_faces = int(sys.argv[1])
vdwfname = sys.argv[2]

# Build (empty) objects
vdw = FFEA_vdw.FFEA_vdw(vdwfname)
vdw.set_num_faces(num_faces)
vdw.write_to_file(vdwfname)
