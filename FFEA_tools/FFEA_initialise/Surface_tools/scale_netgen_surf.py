import sys
import FFEA_surf
if len(sys.argv) != 4:
	sys.exit("Usage: python scale_netgen_surf.py [INPUT .surf] [OUTPUT .surf] [scale]")

# Get args
insurffname = sys.argv[1]
outsurffname = sys.argv[2]
scale = float(sys.argv[3])

# Make and scale a surf
surf = FFEA_surf.FFEA_surf(insurffname)
surf.scale(scale)

# Write to out_fname
surf.write_to_netgen_surf(outsurffname)


