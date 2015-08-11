#!/usr/bin/python

import sys, re, getopt, math
################################################################################
# HELP
################################################################################
HELP = """Data format converter from ffea+trajectory.out to VTK
\t-h or --help : print this help
\t-s : requires argument with path to .ffea file\n
"""
#input = ffea
#			- .top
#			- .node
#			- .surf
#			- .trajectory

################################################################################
# FUNCTIONS
################################################################################
##prepare description of mesh structure for write
# @param inputList (list)	- list of input elements
# @param cellsNum (int)		- number of numberes describing each element
# @param typesNum (string)	- VTK type number
# examples:
#	prepareMeshStructure(elementList1, 4, "10") #mesh
#	prepareMeshStructure(surfList, 3, "5") #meshSurf
def prepareMeshStructure(inputList, cellsNum, typesNum):
	#cells
	meshStr = "\nCELLS "+str(len(inputList))+" "+str((cellsNum+1)*len(inputList))+"\n"
	for i in inputList:
		meshStr += str(cellsNum)+" "+" ".join(i)+"\n"

	#cell_types
	meshStr +="\nCELL_TYPES "+str(len(inputList))+"\n"
	typesNum +=" "
	for i in range(len(inputList)):
		meshStr += typesNum
	meshStr +="\n"
	return meshStr

################################################################################
## write MESH / MESH_SURF
# @param name (string)		- name of output file
# @param step (int)			- #step
# @param meshStr (string)	- meshStr is output prepareMeshStructure function
# @param inputList1 (list)	- list of nodes
# examples:
#	write mesh 	- writeMesh("mesh",#step, meshStr, nodeList1, nodeList2)
#	write meshSurf - writeMesh("meshSurf",#step, meshStr, nodeList1,[])
def writeMesh(name, step, meshStr, inputList1, inputList2):
	outputfile = open (name+str(step)+".vtk", 'w')
	#header
	outputfile.write("# vtk DataFile Version 3.1\n"+
					"*description-name of input file* step:"+str(step)+" \n"+
					"ASCII\nDATASET UNSTRUCTURED_GRID\n"+
					"POINTS "+str(len(inputList1)+len(inputList2))+" FLOAT\n")
	#points
	for i in inputList1: #print ', '.join(str(x) for x in list_of_ints)
		outputfile.write(" ".join(str(x) for x in i)+"\n")
	for i in inputList2:
		outputfile.write(" ".join(str(x) for x in i)+"\n")

	#CELLS + CELL_TYPES
	outputfile.write(meshStr)
	outputfile.close()

################################################################################
## parse the input .NODE file
# infor about number of surface and interior nodes
# @param inputfile (string)	- name of .node file
def getNode(inputfile):
	inputfile = open(inputfile,'r')

	state = 0
	numSurf = -1
	numIntr = -1
	c = 0
	for line in inputfile.readlines():
		c = c+1
		if state == 0:	#waiting for surface nodes
			cont = re.match(r"num_surface_nodes\s?(\d+)\s*",line)
			if cont != None:
				numSurf = int(cont.group(1))
				state = 1
		if state == 1:
			cont = re.match(r"num_interior_nodes\s?(\d+)\s*",line)
			if cont != None:
				numIntr = int(cont.group(1))
				break
	inputfile.close()
	if numSurf==-1 or numIntr==-1:
		sys.stderr.write(".NODE file missing content!\n")
		exit(1)
	return [numSurf, numIntr]

################################################################################
## parse the input .TOP file 
# info about triangular pyramids which builds the blob
# @param inputfile (string)	- name of .top file
def getTop(inputfile):
	inputfile = open(inputfile,'r')

	state = 0
	elementList1 = [] #surface elements
	elementList2 = [] #interior elements
	c=0
	for line in inputfile.readlines():
		c =c+1
		if state == 0:	#waiting for surface elements
			cont = re.match(r"surface elements:.*",line)
			if cont != None:
				state = 1
		elif state == 1: #list of surface nodes
			cont = re.match(r"(\d+) (\d+) (\d+) (\d+).*",line)
			if cont == None:
				state = 2
			else:
				elementList1.append([cont.group(1), cont.group(2), cont.group(3), cont.group(4)])
		elif state == 2: #list of interior elements
			cont = re.match(r"(\d+) (\d+) (\d+) (\d+).*",line)
			if cont == None:
				sys.stderr.write("unknown file (.TOP) content! - line "+str(c)+"\n")
				exit(1)
			else:
				elementList2.append([cont.group(1), cont.group(2), cont.group(3), cont.group(4)])
	inputfile.close()
	if elementList1==[] and elementList2==[]:
		sys.stderr.write(".TOP file missing content!\n")
		exit(1)
	return [elementList1, elementList2]

################################################################################
## parse the input .SURF file 
# info about surface triangles
# @param inputfile (string)	- name of .surf file
def getSurf(inputfile):
	inputfile = open(inputfile,'r')

	state = 0
	surfList = [] #surface_faces
	c=0
	for line in inputfile.readlines():
		c = c+1
		if state == 0:	#waiting for surface nodes
			cont = re.match(r"faces:.*",line)
			if cont != None:
				state = 1
		elif state == 1: #list of surface nodes
			cont = re.match(r"\d+ (\d+) (\d+) (\d+)\s*",line)
			if cont == None:
				sys.stderr.write("unknown file (.NODE) content! - line "+str(c)+"\n")
				exit(1)
			else:
				surfList.append([cont.group(1), cont.group(2), cont.group(3)])
	inputfile.close()
	if surfList==[]:
		sys.stderr.write(".SURF file missing content!\n")
		exit(1)
	return surfList
################################################################################
## F = sqrt( F_x^2 + F_y^2 + F_z^2)
# @param force (list)	- list of input vectors of force
# @return (list)		- list of scalar values of force
def cntForce(force):
	res = []
	for f in force:
		res.append(math.sqrt(math.pow(float(f[0]),2)+math.pow(float(f[1]),2)+math.pow(float(f[2]),2)))
	return res
		
################################################################################
## add forces information to points
# @param outputfile (string)	- name of output file
# @param data (list)		- list of forces in each point
"""
def writeData(outputfile, data):
	scalars = cntForce(data)
	outputfile = open(outputfile,'a')
	outputfile.write("\nPOINT_DATA "+str(len(scalars))+"\nSCALARS scalar_force float 1\nLOOKUP_TABLE default\n")
	for i in scalars:
		outputfile.write(str(i)+"\n")
	outputfile.write("VECTORS vector_force float\n")
	for i in data:
		outputfile.write(" ".join([i[0],i[1],i[2],'\n']))
	outputfile.close()
"""	
def writeData(outputfile, data):
	outputfile = open(outputfile,'a')
	outputfile.write("\nPOINT_DATA "+str(len(data))+"\nVECTORS vector_force float\n")
	for i in data:
		outputfile.write(" ".join([i[0],i[1],i[2],'\n']))
	scalars = cntForce(data)
	outputfile.write("\nSCALARS scalar_force float\nLOOKUP_TABLE default\n")
	for i in scalars:
		outputfile.write(str(i)+"\n")
	outputfile.close()



################################################################################
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
################################################################################
# GETOPT
################################################################################
try:
	optlist, args= getopt.getopt(sys.argv[1:], 's:h', ['help'])
except getopt.GetoptError, err:
    print str(err)
    sys.exit(2)

ok=False
for opt in optlist:
	if "--help" in opt or "-h" in opt:
		sys.stderr.write(HELP)
		sys.exit(0)
	if "-s" in opt:
		inputfile = (opt[1])
		ok=True

if ok == False:
	sys.stderr.write(HELP)
	sys.exit(2)

################################################################################
## parse .ffea file
#CONTENT
# general:
# 	<trajectory_out_fname = sphere_63_120_two_vdw/sphere_63_120_two_vdw_trajectory.out>
#	<num_blobs = 2>
# per blob:
#	<nodes = sphere_63_120_structure/sphere_63_120.node>
#	<topology = sphere_63_120_structure/sphere_63_120.top>
#	<surface = sphere_63_120_structure/sphere_63_120.surf>
################################################################################
path = '/'.join(inputfile.split('/')[:-1])+'/'

inputfile = open(inputfile,'r')
inputBlob = []
state = 0
c = 0
for line in inputfile.readlines():
	#general info
	if state == 0:	#trajectory file
		cont = re.match(r"\s*<trajectory_out_fname\s?=\s?(.*\.out).*",line)
		if cont != None:
			inputTraj = path+'/'+cont.group(1)
			state = 1
	elif state == 1: #nuber of blobs
		cont = re.match(r"\s*<num_blobs\s?=\s?(\d+).*",line)
		if cont != None:
			blobCnt = int(cont.group(1))
			for i in range(blobCnt):
				inputBlob.append([])
			blobC = 0
			state = 2
	#info for each blob
	elif state == 2:
		cont = re.match(r"\s*<nodes\s?=\s?(.*\.node).*",line)
		if cont != None:
			inputBlob[blobC].append(path+cont.group(1))
			state = 3
	elif state == 3:
		cont = re.match(r"\s*<topology\s?=\s?(.*\.top).*",line)
		if cont != None:
			inputBlob[blobC].append(path+cont.group(1))
			state = 4
	elif state == 4:
		cont = re.match(r"\s*<surface\s?=\s?(.*\.surf).*",line)
		if cont != None:
			inputBlob[blobC].append(path+cont.group(1))
			state = 2
			blobC +=1
if state != 2:
	sys.stderr.write("Noncomplete .ffea file")
	exit(2)



# list of all information about blobs
blobInfo = []
#<blob>
#	<node>
#		<#surf nodes>
#		<#inter nodes>
#	<top>
#		<surface elements>
#		<interior elements>
#	<surf>
#		<surface faces>
meshInfo = [[],[]]
#<mesh>
#	<mesh description #1>
#	<mesh description #2>
#	<...>
#<meshSurf>
#	<surf description #1>
#	<surf description #2>
#	<...>
blobMeshIndex=[]
################################################################################
# parse files with information about blob structure (.node, .top, .surf)
################################################################################
for blob in range(blobCnt):
	blobMeshIndex.append([])
	blobInfo.append([])
	
	for data in range(len(inputBlob[blob])):
		for i in range(blob): #check if the input is not the same as input of previos blob
			if inputBlob[blob][data] == inputBlob[i][data]: #same type of structure
				blobMeshIndex[blob].append(blobMeshIndex[i][data]) #index to list with descrition
				break
		if len(blobMeshIndex[blob]) == data: # another type of blob structure
			if data == 0:
				#parse NODE
				blobMeshIndex[blob].append(getNode(inputBlob[blob][data]))
			elif data == 1:
				#parse TOP
				lst = getTop(inputBlob[blob][data]) #------------------------------------------ INTERIOR ELEMENTS?
				meshInfo[0].append(prepareMeshStructure(lst[0], 4, "10"))
				blobMeshIndex[blob].append(len(meshInfo[0])-1)
			elif data == 2:
				#parse SURF
				meshInfo[1].append(prepareMeshStructure(getSurf(inputBlob[blob][data]), 3, "5"))
				blobMeshIndex[blob].append(len(meshInfo[1])-1)

################################################################################
#-----------#
#	TEST	#
#-----------#
#for i in range(len(blobMeshIndex)):
#
#print blobMeshIndex[0]
#	print " "
#print meshInfo[1][blobMeshIndex[0][1]]
#exit()
#-----------#
################################################################################
# main part of code
# read position of points in each timestep
# write all info into output file in VTK file format
# input:	trajectory.out
################################################################################
# read the input trajectory.out file
inputfile = open(inputTraj,'r')

state = 0
nodeMove = {} #surface_faces

nodeList1 = [] #surface
nodeList2 = [] #interior
force = {}
blobID = -1
step = 0

c=0
for line in inputfile.readlines():
	c = c+1

	if state == 1: #DYNAMIC / STATIC info
		state = 2
	elif state == 2: #list of node moves
#		cont = re.match(r"(-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+)\s*",line)
		cont = re.match(r"(-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+).*\s(-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+)\s*$",line)
		if cont == None: #another blob - next step
			state = 0
		else:
			pointN +=1
			if pointN <= blobMeshIndex[blobID][0][0]:
				nodeMove[blobID][0].append([cont.group(1), cont.group(2), cont.group(3)])
			else:
				nodeMove[blobID][1].append([cont.group(1), cont.group(2), cont.group(3)])
			force[blobID].append([cont.group(4), cont.group(5), cont.group(6)])
	if state == 0:	#waiting for blob info
		cont = re.match(r"Blob (\d+), Conformation \d+, step (\d+).*",line)
		if cont != None:
			state = 1
			stepRead = int(cont.group(2))
			if step != stepRead:								
				#write
				for i in nodeMove:
					writeMesh("out/"+str(i)+"mesh",step, meshInfo[0][blobMeshIndex[i][1]], nodeMove[i][0], nodeMove[i][1]) #write mesh
					writeData("out/"+str(i)+"mesh"+str(step)+".vtk",force[i]) #write data (force)
					writeMesh("out/"+str(i)+"meshSurf",step, meshInfo[1][blobMeshIndex[i][2]], nodeMove[i][0],[]) #write meshSurf
				step = stepRead
			blobID = int(cont.group(1))
			nodeMove[blobID] = [[],[]]
			force[blobID] = []
			pointN =0
inputfile.close()















