import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python first_frame_to_pdb [INPUT WALRUS TRAJECTORY FILE] [OUTPUT PDB FILE]")

walrus_traj = open(sys.argv[1], "r")
pdb_out = open(sys.argv[2], "w")
pdb_out.write("REMARK Created with first_frame_to_pdb.py\n")
pdb_out.write("MODEL      1\n");

# skip asterisk line
walrus_traj.readline()

# read the blob x, step y line
current_frame = walrus_traj.readline()
if len(current_frame) == 0:
	sys.exit("No blob x, stp y line...");
print "Converting 1st frame: " + current_frame

#STATIC etc
walrus_traj.readline()

# get number of nodes
num_nodes = int(walrus_traj.readline())
print "num_nodes = " + str(num_nodes)

eof = 0;
while eof == 0:
	for i in range(num_nodes):
		line = walrus_traj.readline()
		if len(line) == 0:
			print "Reached end of file before end of snapshot " + current_frame
			eof = 1
			break
		a = line.split()

		ridiculous_pos_x = ("%.3f" % (float(a[0]) * 1e9)).rjust(12, " ")
		ridiculous_pos_y = ("%.3f" % (float(a[1]) * 1e9)).rjust(8, " ")
		ridiculous_pos_z = ("%.3f" % (float(a[2]) * 1e9)).rjust(8, " ")

		stupid_number_format = str(i).rjust(7, " ")
		pdb_out.write("ATOM" + stupid_number_format + "  N   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
	eof = 1	
print "done"
