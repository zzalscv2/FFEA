import sys, os
from math import *
import FFEA_script

if len(sys.argv) != 4:
	sys.exit("usage: python " + os.path.basename(sys.argv[0]) + " [FFEA script (.ffea) (must contain single blob)] [Number of blobs] [Volume fraction of system]")
	
# Get args
orig_script_fname = sys.argv[1]
num_blobs = int(sys.argv[2])
vol_frac = float(sys.argv[3])

# Read script
script = FFEA_script.FFEA_script(orig_script_fname)
if script.params.num_blobs != 1:
	sys.exit("Error. Script must contain only one blob which will then be repeated to desired number density.")

# Get appropriate params to calc box dimensions
es_h = script.params.es_h
kappa = script.params.kappa


# Get appropriate params to calculate blob dimensions
node_fname = script.blob[0].conformation[0].nodes
top_fname = script.blob[0].conformation[0].topology
scale = script.blob[0].scale

# Get nodes to calculate centroid and maximum radius
nodes = FFEA_nodes.FFEA_nodes(node_fname)

# Build system to calculate volume
topology = FFEA_topology.FFEA_topology(topology_fname)
blob_volume = topology.get_volume(nodes)
total_blob_volume = blob_volume * num_blobs

# For given volume fraction, calculate required box dimensions
bdv = total_blob_volume / vol_frac
bdx = bdy = bdz = math.pow(bdv, 1.0/3.0)
es_N_x = es_N_y = es_N_z = ceil(bdx * kappa / es_h)
bdx = bdy = bdz = es_N_x * 1.0/kappa * es_h
bdv = bdx * bdy * bdz

# Ask user if rounding is ok. If not, rescale blobs a bit
approx_vol_frac = total_blob_volume / bdv
success = raw_input("For your given parameters, volume fraction is %e. Rescale spheres for exact match (y/n)?:" % (approx_vol_frac))
if success.lower == "y":
	new_blob_volume = vol_frac * bdv / num_blobs
	script.blob[0].scale *= pow(new_blob_volume / blob_volume, 1.0/3.0)

# Now add a load of new blobs to the script arranged randomly throughout the box
for i in range(num_blobs - 1):
	script.blob.append(blob[0])

for i in range(0, num_blobs, 1):
	for j in range(0, i, 1):
		# Make sure they don't overlap!
		script.blob.set_centroid_pos()	
