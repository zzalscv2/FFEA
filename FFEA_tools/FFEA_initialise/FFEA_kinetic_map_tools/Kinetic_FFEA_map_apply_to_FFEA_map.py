import sys, os
import numpy as np
import FFEA_kinetic_map

if len(sys.argv) != 5:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT MAP A] [INPUT MAP B] [OUTPUT AB] [OUTPUT BA]")

# Get args
mapa_fname = sys.argv[1]
mapb_fname = sys.argv[2]
outmapab_fname = sys.argv[3]
outmapba_fname = sys.argv[4]

# Get maps
mapa = FFEA_kinetic_map.FFEA_kinetic_map(mapa_fname)
mapb = FFEA_kinetic_map.FFEA_kinetic_map(mapb_fname)

# Apply maps
mapab = mapa.apply_to_map(mapb)
mapba = mapb.apply_to_map(mapa)

# Write maps
mapab.write_to_file(outmapab_fname)
mapba.write_to_file(outmapba_fname)
