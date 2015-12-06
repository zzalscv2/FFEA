import sys, os
import numpy as np
import FFEA_binding_sites, FFEA_node, FFEA_surface

if len(sys.argv) != 6:
	sys.exit("Usage " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT FFEA .node file] [INPUT FFEA .surf file] [OUTPUT FFEA .bsites file] [INPUT node to pin around] [INPUT pin radius]")

# Get args
nodefname = sys.argv[1]
surffname = sys.argv[2]
bsitesfname = sys.argv[3]
origin = int(sys.argv[4])
radius = float(sys.argv[5])

# Build objects
bsites = FFEA_binding_sites.FFEA_binding_sites(bsitesfname)
node = FFEA_node.FFEA_node(nodefname)
surf = FFEA_surface.FFEA_surface(surffname)

newbsite = FFEA_binding_sites.FFEA_binding_site("")
faces = []

# Build face list by seeing which faces are within the region
ocentroid = node.pos[origin]
for f in surf.face:
	fcentroid = f.calc_centroid(node)
	if np.linalg.norm(fcentroid - ocentroid) < radius:
		
		# Within radius
		faces.append(surf.face.index(f))

# Build binding_site
bsite.set_structure(faces)
bsites.add_site(bsite)
bsites.write_to_file("test.bsites")
