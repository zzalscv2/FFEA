import sys, os
import FFEA_trajectory
from Vectors import *

# Import trajectory
mytraj = FFEA_trajectory.FFEA_trajectory(sys.argv[1], 285)

# Output stuff
outfile = open(sys.argv[2], "w")
outfile.write("Frame\tSeparation Score\tOverlap Score\tOrientation Score\n\n")
# Scale node file by 1e-10
mytraj.blob[0].scale_frames(1e-10)

# Calc all individual centroids
mytraj.calc_centroid(0)

# Define motor domains relative to MT
sub_blob_origin = mytraj.blob[0].frame[0].centroid
mytraj.define_sub_blob_radially(1, 0, sub_blob_origin, 250e-10, 1000e-10)
mytraj.define_sub_blob_radially(2, 0, sub_blob_origin, 250e-10, 1000e-10)

# Define constant vectors

# n2 is approximately parallel to stalks
#n2 = mytraj.define_vector_through_nodes(0, 1, 2468, 2430) + mytraj.define_vector_through_nodes(0, 2, 2468, 2430)
n2 = vector3(0.9438, 0.3302, 0.0)
n2.normalise()

# n3 is through microtubule
#n3 = mytraj.define_vector_normal_to_plane(0, 0, 677, 500, 983) 
n3 = vector3(0.0, 0.0, 1.0)
n3.normalise()

# n1 is ideal vector when motor domains are overlapping. Cross priduct of other 2
#n1 = vec3_cross_prod(n2, n3)
n1 = vector3(0.3302, -0.9438, 0.0)
n1.normalise()

# Beginning analysis of each frame
for i in range(mytraj.num_frames):
	position_score = 0.0
	orientation_score = 0.0

	# Define motor domain separation vector
	b12 = mytraj.calc_sub_blob_centroid(i, 1, 0) - mytraj.calc_sub_blob_centroid(i, 2, 0)
	
	# Calculate separation score and normalise by maximum posible separation (half motor + stalk + gap + stalk + half motor)
	separation_score = 1 - sqrt(pow(vec3_dot_prod(b12, n2), 2) + pow(vec3_dot_prod(b12, n3), 2))/45.78e-9

	# Calculate and normalise position score
	position_score = sqrt(pow(vec3_dot_prod(b12, n2), 2) + pow(vec3_dot_prod(b12, n3), 2))
	if position_score > 10e-9:
		position_score = 0
	else:
		position_score = 1 - position_score/10e-9
	
	# Calculate vectors through motor domains centers
	b1 = mytraj.define_vector_normal_to_plane(i, 1, 515, 504, 135) 
	b1.normalise()
	b2 = mytraj.define_vector_normal_to_plane(i, 2, 515, 504, 135) 
	b2.normalise()

	# Calculate and normalise orientation score
	vecsum = b1 + b2
	vecsum.scale(0.5)
	orientation_score = vec3_dot_prod(vecsum, n1)

	# Output stuff
	outfile.write(str(i) + "\t" + str(separation_score) + "\t\t" + str(position_score) + "\t" + str(orientation_score) + "\n")

outfile.close()
os.system("gedit " + sys.argv[2] + " &")
