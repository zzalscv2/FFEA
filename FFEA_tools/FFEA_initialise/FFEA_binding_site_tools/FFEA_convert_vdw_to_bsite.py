import sys, os
import FFEA_vdw, FFEA_binding_sites

if len(sys.argv) != 3:
	sys.exit("python FFEA_convert_vdw_to_bsite.py [INPUT .vdw fname] [OUTPUT .bsites fname]")

# Get args
vdwfname = sys.argv[1]
bsitesfname = sys.argv[2]

# Build objects
vdw = FFEA_vdw.FFEA_vdw(vdwfname)
bsites = FFEA_binding_sites.FFEA_binding_sites(bsitesfname)

# Get new bsite
bsite = FFEA_binding_sites.FFEA_binding_site()

faces = []
for i in range(len(vdw.vdw_index)):
	if vdw.vdw_index[i] != -1:
		faces.append(i)

bsite.set_type(0)
bsite.set_structure(faces)
bsites.add_site(bsite)
bsites.write_to_file(bsitesfname)

