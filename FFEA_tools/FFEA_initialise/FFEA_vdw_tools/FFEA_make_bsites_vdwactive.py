import sys, os
import FFEA_vdw, FFEA_binding_sites

if len(sys.argv) != 5:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .bsites] [INPUT .vdw] [OUTPUT .vdw] [vdw type]")

# Get args
bsitefname = sys.argv[1]
vdwinfname = sys.argv[2]
vdwoutfname = sys.argv[3]
vdw_type = int(sys.argv[4])

# Build objects
bsites = FFEA_binding_sites.FFEA_binding_sites(bsitefname)
vdw = FFEA_vdw.FFEA_vdw(vdwinfname)

# Set all bsites as vdw active
for s in bsites.bsites:
	for f in s.face_index:
		vdw.vdw_index[f] = vdw_type

# Write new file
vdw.write_to_file(vdwoutfname)
