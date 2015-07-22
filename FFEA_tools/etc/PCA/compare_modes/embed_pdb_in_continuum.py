import os, sys

def read_walrus_traj_frame(file):
	# skip asterisk
	check = file.readline()

	if check == "":
		return -1, None

	# skip Blob... Step... line
	file.readline()

	# read the number of nodes
	num_nodes = int(file.readline())

	frame = []
	for i in xrange(num_nodes):
		line = (file.readline()).split()
		frame.append([float(line[0]), float(line[1]), float(line[2])])

	return num_nodes, frame

def read_pdb_frame(file):
	while True:
		line = file.readline()
		if line == "":
			break
		line = line.replace("-", " -")
		line = line.split()
		if "REMARK" in line[0]:
			continue
		if "MODEL" in line[0]:
			continue
		if "ATOM" in line[0]:
			frame = []
			while True:
				if "ENDMDL" in line[0]:
					break
				if "TER" in line[0]:
					break
				if "END" in line[0]:
					break
				frame.append([float(line[5]), float(line[6]), float(line[7])])
				line = (file.readline())
				line = line.replace("-", " -")
				line = line.split()
			return len(frame), frame
	return -1, None

def invert_matrix(M):
	m_inv = [[0,0,0],[0,0,0],[0,0,0]]
        m_inv[0][0] = M[2][2]*M[1][1] - M[2][1]*M[1][2];
        m_inv[0][1] = M[2][1]*M[0][2] - M[2][2]*M[0][1];
        m_inv[0][2] = M[1][2]*M[0][1] - M[1][1]*M[0][2];
        m_inv[1][0] = M[2][0]*M[1][2] - M[2][2]*M[1][0];
        m_inv[1][1] = M[2][2]*M[0][0] - M[2][0]*M[0][2];
        m_inv[1][2] = M[1][0]*M[0][2] - M[1][2]*M[0][0];
        m_inv[2][0] = M[2][1]*M[1][0] - M[2][0]*M[1][1];
        m_inv[2][1] = M[2][0]*M[0][1] - M[2][1]*M[0][0];
        m_inv[2][2] = M[1][1]*M[0][0] - M[1][0]*M[0][1];

        det = M[0][0] * m_inv[0][0] + M[1][0] * m_inv[0][1] + M[2][0] * m_inv[0][2];

        m_inv[0][0]/=det
	m_inv[0][1]/=det
	m_inv[0][2]/=det;
        m_inv[1][0]/=det;
	m_inv[1][1]/=det;
	m_inv[1][2]/=det;
        m_inv[2][0]/=det;
	m_inv[2][1]/=det;
	m_inv[2][2]/=det;

	return m_inv

def mat_mult(M, v):
	result = [0, 0, 0]
	for i in xrange(3):
		for j in xrange(3):
			result[i] += M[i][j]*v[j]
	return result

def get_barycentric_coords(point, t1, t2, t3, t4):
	T = [[t1[0] - t4[0], t2[0] - t4[0], t3[0] - t4[0]], [t1[1] - t4[1], t2[1] - t4[1], t3[1] - t4[1]], [t1[2] - t4[2], t2[2] - t4[2], t3[2] - t4[2]]]
	
	T_inv = invert_matrix(T)

	point_minus_t4 = [0, 0, 0]
	for i in xrange(3):
		point_minus_t4[i] = point[i] - t4[i]

	result = mat_mult(T_inv, point_minus_t4)

	return result[0], result[1], result[2]

def r2_to_tet_centre(point, t1, t2, t3, t4):
	centroid = [(t1[0] + t2[0] + t3[0] + t4[0])/4, (t1[1] + t2[1] + t3[1] + t4[1])/4, (t1[2] + t2[2] + t3[2] + t4[2])/4]
	dx = point[0] - centroid[0]
	dy = point[1] - centroid[1]
	dz = point[2] - centroid[2]
	return dx * dx + dy * dy + dz * dz

def what_tet_is_this_point_in(point, num_elements, topology, traj_frame):
	for j in xrange(num_elements):
		a, b, c = get_barycentric_coords(point, traj_frame[topology[j][0]], traj_frame[topology[j][1]], traj_frame[topology[j][2]], traj_frame[topology[j][3]])
		if a < 0 or a > 1 or b < 0 or b > 1 or c < 0 or c > 1:
			continue
		return [j, a, b, c]

	print "Warning. One or more points are outside of all tets. Picking closest tet..."
	nearest_tet = -1
	shortest_dist2 = float("inf")
	for j in xrange(num_elements):
		dist2 = r2_to_tet_centre(point, traj_frame[topology[j][0]], traj_frame[topology[j][1]], traj_frame[topology[j][2]], traj_frame[topology[j][3]])
		if dist2 < shortest_dist2:
			shortest_dist2 = dist2
			nearest_tet = j

	a, b, c = get_barycentric_coords(point, traj_frame[topology[nearest_tet][0]], traj_frame[topology[nearest_tet][1]], traj_frame[topology[nearest_tet][2]], traj_frame[topology[nearest_tet][3]])
	return [nearest_tet, a, b, c]

	

def get_pdb_to_continuum_map(num_elements, topology, traj_frame, num_pdb_nodes, pdb_frame):
	mapping = []
	for i in xrange(num_pdb_nodes):
		tet_index = what_tet_is_this_point_in(pdb_frame[i], num_elements, topology, traj_frame)
		mapping.append(tet_index)
	return mapping

def map_point(a, b, c, t1, t2, t3, t4):
	new_x = t4[0] + a * (t1[0] - t4[0]) + b * (t2[0] - t4[0]) + c * (t3[0] - t4[0])
	new_y = t4[1] + a * (t1[1] - t4[1]) + b * (t2[1] - t4[1]) + c * (t3[1] - t4[1])
	new_z = t4[2] + a * (t1[2] - t4[2]) + b * (t2[2] - t4[2]) + c * (t3[2] - t4[2])
	return [new_x, new_y, new_z]

def write_pdb_frame(mapping, topology, traj_frame):
	i = 0
	for point_map in mapping:
		tet = topology[point_map[0]]
		mapped_point = map_point(point_map[1], point_map[2], point_map[3], traj_frame[tet[0]], traj_frame[tet[1]], traj_frame[tet[2]], traj_frame[tet[3]])

		ridiculous_pos_x = ("%.3f" % (mapped_point[0])).rjust(12, " ")
		ridiculous_pos_y = ("%.3f" % (mapped_point[1])).rjust(8, " ")
		ridiculous_pos_z = ("%.3f" % (mapped_point[2])).rjust(8, " ")

		stupid_number_format = str(i).rjust(7, " ")
		pdb_out.write("ATOM" + stupid_number_format + "  N   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
		i += 1

	pdb_out.write("TER\nENDMDL\n")

if len(sys.argv) != 5:
	sys.exit("Usage: python embed_pdb_in_continuum.py [CONTINUUM WALRUS TRAJ FILE] [CONTINUUM TOPOLOGY FILE] [INPUT PDB] [OUTPUT PDB]")

print "Reading topology file " + sys.argv[2] + "..."
topfile = open(sys.argv[2])
line = topfile.readline()
if "walrus topology file" not in line:
	sys.exit("Error: Provided topology file is missing 'walrus topology file' as first line...")

line = (topfile.readline()).split()
num_elements = int(line[1])
print "num_elements =", num_elements

line = (topfile.readline()).split()
num_surface_elements = int(line[1])
print "num_surface_elements =", num_surface_elements

line = (topfile.readline()).split()
num_interior_elements = int(line[1])
print "num_interior_elements =", num_interior_elements

line = topfile.readline()
if "surface elements:" not in line:
	sys.exit("Error: Provided topology file is missing 'surface_elements:' line...")

topology = []
while True:
	line = (topfile.readline()).split()
	if "interior" in line[0]:
		break
	topology.append([int(x) for x in line])

for i in xrange(num_interior_elements):
	line = (topfile.readline()).split()
	topology.append([int(x) for x in line])

topfile.close()

trajfile = open(sys.argv[1], "r")
pdbfile = open(sys.argv[3], "r")

print "Reading first frame of walrus traj file", sys.argv[1]
num_nodes, traj_frame = read_walrus_traj_frame(trajfile)
print "num_nodes =", num_nodes

print "Reading first frame of pdb file", sys.argv[3]
num_pdb_nodes, pdb_frame = read_pdb_frame(pdbfile)
pdbfile.close()

print "Generating pdb to continuum mapping..."
mapping = get_pdb_to_continuum_map(num_elements, topology, traj_frame, num_pdb_nodes, pdb_frame)

print "Writing new pdb file..."
frame = 0
pdb_out = open(sys.argv[4], "w")
print "Mapping frame", frame
write_pdb_frame(mapping, topology, traj_frame)
frame += 1
while True:
	num_nodes, traj_frame = read_walrus_traj_frame(trajfile)
	if num_nodes == -1:
		break
	print "Mapping frame", frame
	write_pdb_frame(mapping, topology, traj_frame)
	frame += 1

trajfile.close()
pdb_out.close()

print "Done."
