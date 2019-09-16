import sys, os
import numpy as np
import copy

import MDAnalysis as mda
import FFEA_node, FFEA_trajectory
from scipy.linalg import expm

import icp


def eulerRotate(v, xyz):

	''' Rotate vector v (or array of vectors) by the euler angles xyz '''
	# https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
	for theta, axis in zip(xyz, np.eye(3)):
		v = np.dot(np.array(v), expm(np.cross(np.eye(3), axis*-theta)))
	return v

def printProgressBar(start_text, num_iterations, num_done):

	"""
	Badly print an ascii progress bar.
	"""
	try:
		rows, cols = os.popen('stty size', 'r').read().split()
		cols = int(cols)
	except ValueError:
		cols = 80, 80

	ratio_done = float(num_done)/float(num_iterations)
	str_done = str(num_done)+"/"+str(num_iterations)+" "
	textstr = start_text+" ["
	bar_length= cols-len(textstr.decode('utf-8')+"]")-len(str_done.decode('utf-8'))
	blank_chars = "▒"*int((bar_length*(1-ratio_done)))
	filled_chars = "█"*int((bar_length*(ratio_done)))
	spacer = ""
	str_to_write = "\r"+textstr+filled_chars+blank_chars+"] "+str_done

	if len(str_to_write.decode('utf-8')) > cols:
		blank_char_num_new = int((bar_length*(1-ratio_done)))-(len(str_to_write.decode('utf-8'))-cols)
		blank_chars = "▒"*blank_char_num_new
		str_to_write = "\r"+textstr+filled_chars+blank_chars+"] "+str_done

	if len(str_to_write.decode('utf-8')) < cols:
		spacer = " "*(cols-len(str_to_write.decode('utf-8')))

	if len(str_to_write.decode('utf-8')) != cols:
		return

	sys.stdout.write(str_to_write+spacer)

def calculateCentroid(pos):
	return np.mean(pos, axis=0)

def alignTranslate(pointClouds):

	cent = []
	for pCloud in pointClouds:
		cent.append(calculateCentroid(pCloud))

	trans = cent[1] - cent[0]

	return trans

def alignRotate(pointClouds, numCandidates):

	pcOut = {}

	for i in range(numCandidates):
#		printProgressBar("Calculating optimal alignment", numCandidates, i)
		tempPC = copy.copy(pointClouds[0])
		XYZ = [np.random.random()*2*np.pi, np.random.random()*2*np.pi, np.random.random()*2*np.pi]

		print("Calculating optimal alignment for candidate %d: (%.2f,%.2f,%.2f)" % (i, XYZ[0], XYZ[1], XYZ[2]))

		tempPC = eulerRotate(tempPC, XYZ)
		T, distances = icp.icp(tempPC, pointClouds[1], max_iterations=1000, tolerance=0.01)

		pcOut[np.average(distances)] = [XYZ, T]
	
	return min(pcOut.keys()), pcOut[min(pcOut.keys())][0], pcOut[min(pcOut.keys())][1]

def applyTransformation(pos, T):
    """
    Apply a 4x4 transformation matrix to our 2-D array of 3-D points. Returns
    the newly translated array.
    """
    transArray = np.zeros([ len(pos), 3 ])
    for i in range(len(pos)):
        node_1 = [pos[i][0], pos[i][1], pos[i][2], 1]
        transArray[i] = np.dot(T, node_1)[:3]
    return transArray

# Get args (manually)
if len(sys.argv) < 3:
	sys.exit("Usage: python " + os.path.basename(sys.argv[0]) + " [Structure 1 (.node)] [Structure 2 (.pdb)] Optional [numCandidates]")

iFnames = [sys.argv[1], sys.argv[2]]

try:
	numCandidates = int(sys.argv[3])
except:
	numCandidates = 1

obj = []
pc = []

# Check consistency and read in
for f in iFnames:

	# Exists
	if(not os.path.exists(f)):
		sys.exit("Error: " + f + " could not be found :( ")

	# Extension
	base, ext = os.path.splitext(f)
	
	if(ext == ".node"):
		obj.append(FFEA_node.FFEA_node())

		try:
			obj[-1].load(f)
			pc.append(obj[-1].pos)
		except:
			sys.exit("Error: " + f + " could not be read from file :( ")

	elif(ext == ".pdb"):
		try:
			obj.append(mda.Universe(f))
			C = obj[-1].select_atoms("protein")
			pc.append(C.positions)
		except:
			sys.exit("Error: " + f + " could not be read from file :( ")


#
# Align (1 onto 2)
#

# Translation
trans = alignTranslate(pc)
pc[0] += trans

# Rotation
rmsd, rotAng, rotMat = alignRotate(pc, numCandidates)

print("Optimal alignment found:")
print("Translation: %.2f, %.2f %.2f" % (trans[0], trans[1], trans[2]))
print("Rotation: %.2f, %.2f %.2f" % (rotAng[0], rotAng[1], rotAng[2]))
#print("Rotation matrix: \n" +str(rotMat))
print("RMSD: " + str(rmsd))

# Apply candidate rotation and transformation
pc[0] = eulerRotate(pc[0], rotAng)
pc[0] = applyTransformation(pc[0], rotMat)

# Write
iFname = iFnames[0]
obj = obj[0]

# Extension
base, ext = os.path.splitext(iFname)
oFname = base + "_Aligned" + ext

if(ext == ".node"):

	obj.pos = pc[0]
	obj.write_to_file(oFname)

elif(ext == ".pdb"):
	sys.exit("Currently unsupported output :(")


