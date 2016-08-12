import os, sys

if len(sys.argv) != 9:
	sys.exit("Usage: python create_homo_material_params_file.py [NUM_ELEMENTS] [DENSITY] [SHEAR VISC] [BULK VISC] [SHEAR MOD] [BULK MOD] [DIELECTRIC] [OUTPUT FILE]")

num_elem = int(sys.argv[1])
density = float(sys.argv[2])
shear_visc = float(sys.argv[3])
bulk_visc = float(sys.argv[4])
shear_modulus = float(sys.argv[5])
bulk_modulus = float(sys.argv[6])
dielec = float(sys.argv[7])

outfile = open(sys.argv[8], "w")
outfile.write("walrus material params file\n")
num_elements = 0
s = ''
for i in xrange(num_elem):
	s += str(density) + " " + str(shear_visc) + " " + str(bulk_visc) + " " + str(shear_modulus) + " " + str(bulk_modulus) + " " + str(dielec) + "\n"

outfile.write("num_elements " + str(num_elements) + "\n")
outfile.write(s)
outfile.close()
