import sys, os
import FFEA_pdb

if len(sys.argv) != 3:
	sys.exit("Usage: python get_first_frame.py [INPUT .pdb] [OUTPUT .pdb]")

pdb = FFEA_pdb.FFEA_pdb(sys.argv[1], num_frames_to_read = 1)
pdb.write_to_file(sys.argv[2])

