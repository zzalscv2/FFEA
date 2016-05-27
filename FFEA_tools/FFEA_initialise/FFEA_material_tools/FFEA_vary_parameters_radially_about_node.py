import FFEA_material, FFEA_node, FFEA_topology
import sys, os
import numpy as np

if len(sys.argv) != 13:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .mat file] [OUTPUT .mat file] [INPUT .node file] [INPUT .top file] [new parameter radius] [central node] [density] [shear_visc] [bulk_visc] [shear_mod] [bulk_mod] [dielec]\n")

# Get args
inmatfname = sys.argv[1]
outmatfname = sys.argv[2]
innodefname = sys.argv[3]
intopfname = sys.argv[4]
radius = float(sys.argv[5])
node_index = int(sys.argv[6])
d = float(sys.argv[7])
sv = float(sys.argv[8])
bv = float(sys.argv[9])
sm = float(sys.argv[10])
bm = float(sys.argv[11])
di = float(sys.argv[12])

# Build objects
mat = FFEA_material.FFEA_material(inmatfname)
node = FFEA_node.FFEA_node(innodefname)
top = FFEA_topology.FFEA_topology(intopfname)

# Run sweep and set parameters
central_node = node.pos[node_index]
num_changed = 0
for el in top.element:
	elindex = top.element.index(el)
	if np.linalg.norm(el.calc_centroid(node) - central_node) < radius:
		mat.set_params(elindex, d, sv, bv, sm, bm, di)
		num_changed += 1

print "Changed " + str(num_changed) + " elements, leaving " + str(top.num_elements - num_changed) + " elements as they were."
mat.write_to_file(outmatfname)
	
