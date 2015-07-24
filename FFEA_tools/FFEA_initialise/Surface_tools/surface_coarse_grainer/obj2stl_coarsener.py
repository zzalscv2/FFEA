import sys, os

argc = len(sys.argv)
if argc != 5:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT input_fname (.surf)] [OUTPUT output_fname (.stl)] [COARSENING TYPE FLAG] [COARSENING LEVEL]")

in_fname = sys.argv[1]
out_fname = sys.argv[2]

i = 3
while(i < argc):
	flag = sys.argv[i]
	i += 1
	if(flag == "-h"):
		break
	level = sys.argv[i]
	print level, flag
	i += 1
	if flag == "-l":
		print "Coarsening by length..."
	elif flag == "-n":
		print "Coarsening by number of edges..."

gts_fname = "gts_fname.gts"
stl_fname = "stl_fname.stl"
coarse_gts_fname = "coarse_gts_fname.gts"

# First obj2gts
os.system("python /localhome/py09bh/FFEA/development/FFEA_git/FFEA_tools/convert_scripts/surf2gts.py " + in_fname + " " + gts_fname)

# gts2stl stl2gts (no idea why)
os.system("gts2stl < " + gts_fname + " > " + stl_fname)
os.system("stl2gts < " + stl_fname + " > " + gts_fname)

# Then gts coarsen
os.system("~/Downloads/gts-0.7.6/examples/coarsen " + flag + " " + level + " -v < " + gts_fname + " > " + coarse_gts_fname)

if(flag == "-h"):
	sys.exit()

# Finally gts2stl
os.system("~/Downloads/gts-0.7.6/tools/gts2stl < " + coarse_gts_fname + " > " + out_fname)

# Get rid of temp files
#os.system("rm " + gts_fname + " " + coarse_gts_fname + " " + stl_fname)

