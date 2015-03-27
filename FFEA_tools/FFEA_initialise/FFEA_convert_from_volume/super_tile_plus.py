import math
import os
import sys
import random

def bastard(y, z):
	b = 2 * (y%2) - 1
	c = 2 * (z%2) - 1
	return .5 - ((b * c) + 1)/2


if len(sys.argv) != 17:
	sys.exit("Usage: python super_tile_plus.py [BLOB FILES STEM NAME] [N_x] [N_y] [N_z] [es_h] [kappa] [scale] [margin_x] [margin_y] [margin_z] [num_blobs (0 for max)] [vdw_eps] [bulk_mod] [shear_mod] [dt] [stagger ('y' or 'n')]")

# path and template name of blob we want to tile
blob = sys.argv[1]

# box dimensions
N_x = int(sys.argv[2])
N_y = int(sys.argv[3])
N_z = int(sys.argv[4])
es_h = float(sys.argv[5])
kappa = float(sys.argv[6])
scale = float(sys.argv[7])
margin_x = float(sys.argv[8])
margin_y = float(sys.argv[9])
margin_z = float(sys.argv[10])

num_blobs = int(sys.argv[11])

vdw_eps = float(sys.argv[12])
bulk_mod = float(sys.argv[13])
shear_mod = float(sys.argv[14])

dt = float(sys.argv[15])

stagger = sys.argv[16]
if stagger != 'y' and stagger !='n':
	sys.exit("Yaw bruv, put in the roight fockin' staggah value, you gettin' me?")

#print "Calculating dimensions of " + blob + ".node"
os.system("python ~/walrus_2013_March_7th/walrus_2013_March_7th/utils/get_dimensions.py " + blob + ".node | grep 'size' > __super_tile_dimensions_output__")

dimensions = open("__super_tile_dimensions_output__", "r")
sizes = dimensions.readlines()
dimensions.close()

size_x = float((sizes[0].split())[2]) * scale
size_y = float((sizes[1].split())[2]) * scale
size_z = float((sizes[2].split())[2]) * scale

max_size = max([size_x, size_y, size_z])
size_x = size_y = size_z = max_size


#print "size_x = ", size_x
#print "size_y = ", size_y
#print "size_z = ", size_z

fpath, fname = os.path.split(blob)

print "<param>"
print "	<restart = 0>"
print ""
print "	<dt = " + str(dt) + ">"
print "	<kT = 4e-21>"
print "	<check = 1000>"
print "	<num_steps = 1e9>"
print "	<rng_seed = time>"
print ""
print "	<trajectory_out_fname = " + fname + "_trajectory.out>"
print "	<measurement_out_fname = " + fname + "_measurement.out>"
print ""
print "	<epsilon = 1e-14>"
print "	<max_iterations_cg = 1000>"
print ""
print "	<kappa = " + str(kappa) + ">"
print ""
print "	<epsilon_0 = 1>"
print "	<dielec_ext = 1>"
print ""
print "	<es_update = 1>"
print "	<es_N_x = " + str(N_x) + ">"
print "	<es_N_y = " + str(N_y) + ">"
print "	<es_N_z = " + str(N_z) + ">"
print "	<es_h = " + str(es_h) + ">"
print ""
print "	<calc_vdw = 0>"
print "	<vdw_r_eq = 1e-9>"
print "	<vdw_eps = " + str(vdw_eps) + ">"
print "	<sticky_wall_xz = 1>"
print "	<calc_es = 0>"
print ""
print "	<wall_x_1 = PBC>"
print "	<wall_x_2 = PBC>"
print "	<wall_y_1 = HARD>"
print "	<wall_y_2 = HARD>"
print "	<wall_z_1 = PBC>"
print "	<wall_z_2 = PBC>"
print ""
print "	<do_stokes = 1>"


cx = N_x * es_h * (1.0/kappa) - 2 * margin_x
cy = N_y * es_h * (1.0/kappa) - 2 * margin_y
cz = N_z * es_h * (1.0/kappa) - 2 * margin_z

num_x = int(math.floor(cx/size_x))
num_y = int(math.floor(cy/size_y))
num_z = int(math.floor(cz/size_z))

if num_blobs == 0:
	num_blobs = num_x * num_y * num_z

print "	<num_blobs = " + str(num_blobs) + ">"
print ""
print "</param>"
print "<system>"

# space out the blobs more sensibly (use all available space)
size_x_ample = cx/num_x
size_y_ample = cy/num_y
size_z_ample = cz/num_z

#print size_x_ample, size_y_ample, size_z_ample

if stagger == 'y':
	i = 0
	for z in range(num_z):
		for y in range(num_y):
			for x in range(num_x):

				print "	<blob>"
				print "		<nodes = " + blob + ".node>"
				print "		<topology = " + blob + ".top>"
				print "		<surface = " + blob + ".surf>"
				print "		<solver = CG>"
				print "		<centroid_pos=(" + str((x + .25 + bastard(y,z)/2.0) * size_x_ample + margin_x) + "," + str((y + .25 + (.5 * (z%2))) * size_y_ample + margin_y) + "," + str((z + .25) * size_z_ample + margin_z) + ")>"
				print "		<velocity=(0,0,0)>"
				print "		<rho = 1.5e3>"
				print "		<shear_visc = 1e-3>"
				print "		<bulk_visc = 1e-3>"
				print "		<shear_mod = " + str(shear_mod) + ">"
				print "		<bulk_mod = " + str(bulk_mod) + ">"
				print "		<scale = " + str(scale) + ">"
				print "	</blob>"	
				i += 1
				if i == num_blobs:
					print "</system>"
					sys.exit(1)
	print "</system>"
else:
	i = 0
	for z in range(num_z):
		for y in range(num_y):
			for x in range(num_x):

				print "	<blob>"
				print "		<nodes = " + blob + ".node>"
				print "		<topology = " + blob + ".top>"
				print "		<surface = " + blob + ".surf>"
				print "		<solver = CG>"
				print "		<centroid_pos=(" + str((x + .5) * size_x_ample + margin_x) + "," + str((y + .5) * size_y_ample + margin_y) + "," + str((z + .5) * size_z_ample + margin_z) + ")>"
				print "		<velocity=(0,0,0)>"
				print "		<rho = 1.5e3>"
				print "		<shear_visc = 1e-3>"
				print "		<bulk_visc = 1e-3>"
				print "		<shear_mod = " + str(shear_mod) + ">"
				print "		<bulk_mod = " + str(bulk_mod) + ">"
				print "		<scale = " + str(scale) + ">"
				print "	</blob>"	
				i += 1
				if i == num_blobs:
					print "</system>"
					sys.exit(1)
	print "</system>"
