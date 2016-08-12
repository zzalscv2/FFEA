import sys
import FFEA_surface, FFEA_node
import FFEA_vdw

if len(sys.argv) != 4:
	sys.exit("Usage: python get_vdw_area.py [INPUT .vdw file] [INPUT .surf file] [INPUT .node file]")

# Get args
vdwfname = sys.argv[1]
surffname = sys.argv[2]
nodefname = sys.argv[3]

# Open files
vdw = FFEA_vdw.FFEA_vdw(vdwfname)
surf = FFEA_surface.FFEA_surface(surffname)
node = FFEA_node.FFEA_node(nodefname)

areas = vdw.calc_active_areas(surf, node)
for i in range(-1,6,1):

	print "VdW index " + str(i) + ": Area = " + str(areas[i + 1]) 
