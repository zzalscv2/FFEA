import sys, os
import numpy as np
import FFEA_binding_sites, FFEA_node, FFEA_surface

if len(sys.argv) != 8:
	sys.exit("Usage " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT FFEA .node file] [INPUT FFEA .surf file] [OUTPUT FFEA .bsites file] [INPUT node to pin around] [INPUT pin radius] [INPUT face angle limit (degrees)] [INPUT site_type]")

# Get args
nodefname = sys.argv[1]
surffname = sys.argv[2]
bsitesfname = sys.argv[3]
origin = int(sys.argv[4])
radius = float(sys.argv[5])
anglelim = (np.pi / 180.0) * float(sys.argv[6])
site_type = int(sys.argv[7])
if site_type < 0:
	sys.exit("Error. Binding site type cannot be less than zero.\n")

# Build objects
bsites = FFEA_binding_sites.FFEA_binding_sites(bsitesfname)
node = FFEA_node.FFEA_node(nodefname)
surf = FFEA_surface.FFEA_surface(surffname)

newbsite = FFEA_binding_sites.FFEA_binding_site()
faces = []
face_indices = []

# Build face list by seeing which faces are within the region
ocentroid = node.pos[origin]
site_centroid = np.array([0.0,0.0,0.0])
for f in surf.face:
	fcentroid = f.calc_centroid(node)
	if np.linalg.norm(fcentroid - ocentroid) < radius:
		
		# Within radius
		faces.append(f)
		face_indices.append(surf.face.index(f))

		site_centroid += fcentroid

# Finish site centroid
site_centroid *= 1.0 / len(faces)

# Find out which face the centroid is in
distance = float("inf")
centindex = None
for f in faces:
	
	if origin not in f.n:
		continue

	fcentroid = f.calc_centroid(node)
	if np.linalg.norm(fcentroid - site_centroid) < distance:
		distance = np.linalg.norm(fcentroid - site_centroid)
		centface = f

# Now check the angle between normals of these faces
centnorm = centface.get_normal(node)
for f in faces:
	if f == centface:
		continue
	
	fnorm = f.get_normal(node)

	# Get rid if angle too big
	if np.arccos(np.dot(centnorm, fnorm)) > anglelim:
		face_indices.remove(surf.face.index(f))

# Build binding_site
newbsite.set_type(site_type)
newbsite.set_structure(face_indices)
bsites.add_site(newbsite)
bsites.write_to_file(bsitesfname)
