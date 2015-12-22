import sys, os
import FFEA_binding_sites, FFEA_node, FFEA_surface
import numpy as np


if len(sys.argv) != 6:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .node file] [INPUT .surf file] [OUTPUT .bsites file] [surface normal vector (x, -x, y, -y, z or -z)] [num_sites_in_row]")

# Get args
external_script = os.path.dirname(os.path.abspath(sys.argv[0])) + "/FFEA_add_binding_site_built_radially_about_node.py"
nodefname = sys.argv[1]
surffname = sys.argv[2]
bsitesfname = sys.argv[3]
norm = sys.argv[4]
num_sites_in_row = int(sys.argv[5])

# Make objects
node = FFEA_node.FFEA_node(nodefname)
surf = FFEA_surface.FFEA_surface(surffname)
bsites = FFEA_binding_sites.FFEA_binding_sites(bsitesfname)
bsites.reset()
bsites.write_to_file(bsitesfname)

# Find object dimensions and scan area
dimens, limits = surf.calculate_structure_dimensions(node)
if norm == "x":
	norm = np.array([1,0,0])
	scandirs = [1,2]
elif norm == "-x":
	norm = np.array([-1,0,0])
	scandirs = [1,2]
elif norm == "y":
	norm = np.array([0,1,0])
	scandirs = [0,2]
elif norm == "-y":
	norm = np.array([0,-1,0])
	scandirs = [0,2]
elif norm == "z":
	norm = np.array([0,0,1])
	scandirs = [0,1]
elif norm == "-z":
	norm = np.array([0,0,-1])
	scandirs = [0,1]

# Get appropriate faces only
faces = []
pos3 = 0
for f in surf.face:
	fnorm = f.get_normal(node)
	if np.arccos(np.dot(norm, fnorm)) * (180 / np.pi) < 10:
		faces.append(f)
		pos3 += f.calc_centroid(node)[1]

pos3 *= 1.0 / len(faces)

# Build binding sites on appropriate surface
num_steps = num_sites_in_row
step = (limits[scandirs[0]][1] - limits[scandirs[0]][0]) / num_steps
radius = step / 3.0

j = 0
for pos1 in [limits[scandirs[0]][0] + (i + 0.5) * step for i in range(num_steps)]:
	sys.stdout.write("\r%5.2f%% completed" % (j * 100.0 / num_steps))
	sys.stdout.flush()
	j += 1
	for pos2 in [limits[scandirs[1]][0] + (i + 0.5) * step for i in range(num_steps)]:
		
		# Find closest node to this position
		pos = np.array([pos1, pos3, pos2])
		min_distance = float("inf")		
		for f in faces:
			for n in f.n:
				distance = np.linalg.norm(node.pos[n] - pos)
				if distance < min_distance:
					min_distance = distance
					central_node = n

		# Add binding site centered at this node
		#os.system("python " + external_script + " " + nodefname + " " + surffname + " " + bsitesfname + " " + str(central_node) + " " + str(radius) + " " + "10" + " " + "0")
		newbsite = FFEA_binding_sites.FFEA_binding_site()
		bsite_faces = []
		bsite_face_indices = []
		for f in faces:
			if np.linalg.norm(f.calc_centroid(node) - pos) < radius:
				bsite_faces.append(f)
				bsite_face_indices.append(surf.face.index(f))
		newbsite.set_type(0)
		newbsite.set_structure(bsite_face_indices)
		bsites.add_site(newbsite)

sys.stdout.write("\r100.00% completed\n")
bsites.write_to_file(bsitesfname)		
