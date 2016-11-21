import sys, os
import FFEA_node, FFEA_pdb

if (len(sys.argv) != 4):
	sys.exit("Usage: python FFEA_map_PDB_to_FFEA.py [INPUT .node] [INPUT .pdb] [INPUT PDB-Node scale]")
	
# Get args
innode = sys.argv[1]
inpdb = sys.argv[2]
scale = float(sys.argv[3])

# Build structures
node = FFEA_node.FFEA_node(innode)
pdb = FFEA_pdb.FFEA_pdb(inpdb)

# Assuming these are not deformed (volume overlap maximisation may be possible, but this is easier for starters)
# Move both to origin
node.set_pos([0.0,0.0,0.0])
pdb.blob[0].frame[0].set_pos([0.0,0.0,0.0])

basepdb = "temp.pdb"
targetnode = "temp.node"
mapfname = "pdbtonode.map"
sparsemapfname = "pdbtonode_sparse.map"
node.write_to_file(targetnode)
pdb.write_to_file(basepdb)

# Make the map
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
os.system(scriptdir + "/make_structure_map -i %s -o %s -m %s -s %f" % (basepdb, targetnode, mapfname, scale))

# Make them sparse
os.system("python " + scriptdir + "/FFEA_convert_kinetic_map_to_sparse.py " + mapfname + " " + sparsemapfname)

os.system("rm %s %s" % (basepdb, targetnode))