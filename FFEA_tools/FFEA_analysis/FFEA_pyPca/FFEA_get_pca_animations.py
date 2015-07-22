import sys, os
import FFEA_traj

if len(sys.argv) != 5:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .pcz file] [INPUT reference topology (_frame0.pdb)] [num_animations] [ffea scale]")

inpcz = sys.argv[1]
inref = sys.argv[2]
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
num_anim = int(sys.argv[3])
anim_basename = inpcz.split(".")[0]
ffea_scale = float(sys.argv[4])

# Get number of nodes for FFEA_stuff
for i in range(num_anim):
	anim_outfname = anim_basename + "_anim" + str(i) + ".pdb"
	anim_outfname_ffea = anim_basename + "_anim" + str(i) + ".out"
	os.system("pyPczdump --input " + inpcz + " --anim " + str(i) + " --pdb " + inref + " -o " + anim_outfname)
	os.system("python " + script_dir + "/../../FFEA_initialise/PDB_tools/PDB_convert_to_FFEA_trajectory/PDB_convert_to_FFEA_trajectory.py " + anim_outfname + " " + anim_outfname_ffea + " " + str(ffea_scale))
