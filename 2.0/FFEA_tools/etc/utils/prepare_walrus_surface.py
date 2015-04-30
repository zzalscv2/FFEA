import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python prepare_walrus_surface.py [OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE] [OUTPUT SURFACE FILE]")

path = os.path.dirname(sys.argv[0])
print "path of running script: " + path

node = sys.argv[1]
top = sys.argv[2]
surf = sys.argv[3]

print "Processing surface..."

# extract the surface
os.system("python " + path + "/extract_surface_from_topology.py " + top + " " + surf + "\n")

# move all surface nodes to the top of the list
os.system("python " + path + "/put_surface_nodes_first.py " + node + " " + top + " " + surf + "\n")

# cycle surface indices so that all normals point out of the mesh
os.system("python " + path + "/make_all_normals_point_outwards.py " + node + " " + top + " " + surf + "\n")

print "Surface preparations complete."
