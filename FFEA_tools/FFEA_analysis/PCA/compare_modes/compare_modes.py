#!/usr/bin/env python
import os, sys

def dot(fname1, fname2):
	eigen1 = open(fname1, "r")
	eigen2 = open(fname2, "r")

	N1 = int((eigen1.readline().split())[2])
	N2 = int((eigen2.readline().split())[2])

	if N1 != N2:
		sys.exit("Vectors are not of equal size (N1 = " + str(N1) + " and N2 = " + str(N2) + " )")

	sum = 0
	for i in xrange(N1):
		e1 = float(eigen1.readline())
		e2 = float(eigen2.readline())
		sum += e1 * e2

	eigen1.close()
	eigen2.close()

	return sum


if len(sys.argv) != 9:
	sys.exit("Usage: ./compare_modes.py [COLON SEPARATED LIST OF 'SMALL MODEL' PDB FNAMES] [COLON SEPARATED LIST OF 'LARGE MODEL' PDB FNAMES] ['LARGE MODEL' TOPOLOGY FILE] [SMALL MODEL REST FRAME] [ALIGN ITERATIONS] [INITIAL YAW] [INITIAL PITCH] [INITIAL ROLL]\nExample: python compare_modes.py s1.pdb:s2.pdb:s3.pdb l1.pdb:l2.pdb:l3.pdb l.top 14 10000 0.0 3.14 0.0\n")

path = os.path.dirname(sys.argv[0])
if path == '':
        path = '.'
print "path of running script: " + path

SMALL_FILES = sys.argv[1].split(':')
LARGE_FILES = sys.argv[2].split(':')

print "'Small model' pdb files given are:"
for f in SMALL_FILES:
	print f

print "'Large model' pdb files given are:"
for f in LARGE_FILES:
	print f

topology = sys.argv[3]
print "'Large model' topology file =", topology

frame_number = sys.argv[4]
num_iter = sys.argv[5]
init_yaw = sys.argv[6]
init_pitch = sys.argv[7]
init_roll = sys.argv[8]

#ENM_FILES=('../atp_anim_007.pdb' '../atp_anim_008.pdb' '../atp_anim_009.pdb' '../atp_anim_010.pdb' '../atp_anim_011.pdb')
#FFEA_FILES=('../emd_1590_05_full_trajectory_anim_1.pdb' '../emd_1590_05_full_trajectory_anim_2.pdb' '../emd_1590_05_full_trajectory_anim_3.pdb' '../emd_1590_05_full_trajectory_anim_4.pdb' '../emd_1590_05_full_trajectory_anim_5.pdb')
#FFEA_TOPOLOGY="../emd_1590_05_full.top"
num_small = len(SMALL_FILES)
num_large = len(LARGE_FILES)

# clear the folders of prior crap
os.system("rm -r __compare_modes_SMALL_TRANSFORMED\n")
os.system("rm -r __compare_modes_LARGE_MODES\n")
os.system("rm -r __compare_modes_LARGE_MAPPED_TO_SMALL\n")
os.system("rm -r __compare_modes_EIGEN\n")

os.system("mkdir __compare_modes_SMALL_TRANSFORMED\n")
os.system("mkdir __compare_modes_LARGE_MODES\n")
os.system("mkdir __compare_modes_LARGE_MAPPED_TO_SMALL\n")
os.system("mkdir __compare_modes_EIGEN\n")

# Get the 14th frame in the ENM pdb (corresponds to the "rest" state)
os.system("python " + path + "/extract_nth_pdb_frame.py " + SMALL_FILES[0] + " " + str(frame_number) + " SMALL_rest.pdb\n")

# Get the translation and rotation required to align the ENM point cloud within the FFEA point cloud
os.system(path + "/../align_point_clouds/align_point_cloud " + LARGE_FILES[0] + " SMALL_rest.pdb aligned_SMALL.pdb align_soln " + num_iter + " " + init_yaw + " " + init_pitch + " " + init_roll + "\n")
align_soln = open("align_soln", "r")
line = align_soln.readline().split()
align_soln.close()
dx = line[0]
dy = line[1]
dz = line[2]
yaw = line[3]
pitch = line[4]
roll = line[5]

# Transform all of the small model pdb modes using this info, storing the transformed pdbs in __compare_modes_SMALL_TRANSFORMED
for small_fname in SMALL_FILES:
	print "Transforming", small_fname
	trans_fname = "__compare_modes_SMALL_TRANSFORMED/" + small_fname + "_TRANSFORMED.pdb"
	os.system("python " + path + "/transform_pdb.py " + small_fname + " " + yaw + " " + pitch + " " + roll + " " + " " + dx + " " + dy + " " + dz + " " + trans_fname + "\n")

# Convert all the FFEA modes into walrus traj files, storing them in __compare_modes_LARGE_MODES
for large_fname in LARGE_FILES:
	print "Converting", large_fname, "to walrus trajectory"
	traj_fname = "__compare_modes_LARGE_MODES/" + large_fname + "_traj.out"
	os.system("python " + path + "/convert_pdb_to_walrus.py " + large_fname + " " + traj_fname + "\n")

# Map the aligned ENM structure onto the frames of the FFEA modes, storing in __compare_modes_LARGE_MAPPED_TO_SMALL
for large_fname in LARGE_FILES:
	traj_fname = "__compare_modes_LARGE_MODES/" + large_fname + "_traj.out"
	mapped_fname = "__compare_modes_LARGE_MAPPED_TO_SMALL/" + large_fname + "_ENM.pdb"
	print "Embedding aligned_ENM.pdb in trajectory", traj_fname, "with topology", topology
	os.system("python " + path + "/embed_pdb_in_continuum.py " + traj_fname + " " + topology + " aligned_SMALL.pdb " + mapped_fname + "\n")

# Get the eigenvectors of the modes in __compare_modes_SMALL_TRANSFORMED
for small_fname in SMALL_FILES:
	trans_fname = "__compare_modes_SMALL_TRANSFORMED/" + small_fname + "_TRANSFORMED.pdb"
	small_eigen_fname = "__compare_modes_EIGEN/" + small_fname + "_eigen"
	print "Calculating eigen vector for", trans_fname
	os.system("python " + path + "/get_eigen_vector.py " + trans_fname + " " + small_eigen_fname + "\n")

# Get the eigenvectors of the modes in __compare_modes_LARGE_MAPPED_TO_SMALL
for large_fname in LARGE_FILES:
	mapped_fname = "__compare_modes_LARGE_MAPPED_TO_SMALL/" + large_fname + "_ENM.pdb"
	large_eigen_fname = "__compare_modes_EIGEN/" + large_fname + "_eigen"
	print "Calculating eigen vector for", mapped_fname
	os.system("python " + path + "/get_eigen_vector.py " + mapped_fname + " " + large_eigen_fname + "\n")

# Dot the ENM and FFEA eigenvectors
dot_prods = [0.0 for i in xrange(num_small * num_large)]
i = 0
for large_fname in LARGE_FILES:
	large_eigen_fname = "__compare_modes_EIGEN/" + large_fname + "_eigen"
	for small_fname in SMALL_FILES:
		small_eigen_fname = "__compare_modes_EIGEN/" + small_fname + "_eigen"
		dot_prods[i] = dot(large_eigen_fname, small_eigen_fname)
		i += 1

dot_prod_out = open("__compare_modes_DOT_PRODUCT_OUTPUT", "w")
for i in xrange(num_small * num_large):
	dot_prod_out.write(str(dot_prods[i]) + "\n")
dot_prod_out.close()

os.system("python " + path + "/make_splot.py __compare_modes_DOT_PRODUCT_OUTPUT " + str(num_small) + " " + str(num_large) + " > __compare_modes_DOT_PRODUCT_OUTPUT_splot\n")
