import sys
import FFEA_surf, FFEA_vdw, FFEA_node

if len(sys.argv) != 4:
	sys.exit("Usage: python get_vdw_area.py [INPUT .vdw file] [INPUT .surf file] [INPUT .node file]")

# Get args
vdwfname = sys.argv[1]
surffname = sys.argv[2]
nodefname = sys.argv[3]

# Open files
vdw = FFEA_vdw.FFEA_vdw(vdwfname)
surf = FFEA_surf.FFEA_surf(surffname)
node = FFEA_node.FFEA_node(nodefname)
surf.load_FFEA_nodes(node)

for i in range(-1,6,1):
	area = FFEA_vdw.get_area_of_index(i, vdw, surf, node)
	print "VdW index " + str(i) + ": Area = " + str(area) 
