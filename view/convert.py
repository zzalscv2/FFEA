#!/usr/bin/python

# first version of script
# only special type of input

import sys, re
################################################################################
# FUNCTIONS
################################################################################
#prepare description of mesh structure for write

#prepareMeshStructure(elementList1, 4, "10") #mesh
#prepareMeshStructure(surfList, 3, "5") #meshSurf

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
#write MESH / MESH_SURF

#write mesh 	- writeMesh("mesh",#step, meshStr, nodeList1, nodeList2)
#write meshSurf - writeMesh("meshSurf",#step, meshStr, nodeList1,[])

# @param name (string)		- name of output file
# @param step (int)			- #step
# @param meshStr (string)	- meshStr is output prepareMeshStructure function
# @param inputList1 (list)	- list of nodes
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
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#fce - radek, pocet cisel, int/float - seznam cisel

################################################################################
## input:	.node, .top
## output:	mesh
################################################################################
# read the input NODE file
inputfile = open("in/input.node",'r')

state = 0
nodeList1 = [] #surface
nodeList2 = [] #interior
c = 0
for line in inputfile.readlines():
	c = c+1
	if state == 0:	#waiting for surface nodes
		cont = re.match(r"surface nodes:.*",line)
		if cont != None:
			state = 1

	elif state == 1: #list of surface nodes
		cont = re.match(r"(-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+) (-?\d+.\d+e[-\+]\d+)\s*",line)
		if cont == None:
			state = 2
		else:
			nodeList1.append([float(cont.group(1)), float(cont.group(2)), float(cont.group(3))])



XXX="""
	elif state == 1: #list of surface nodes
		cont = re.match(r"(-?\d+.\d+) (-?\d+.\d+) (-?\d+.\d+)\s*",line)
		if cont == None:
			state = 2
		else:
			nodeList1.append([float(cont.group(1)), float(cont.group(2)), float(cont.group(3))])
	elif state == 2: #list of interior nodes
		cont = re.match(r"(-?\d+.\d+) (-?\d+.\d+) (-?\d+.\d+)\s*",line)
		if cont == None:
			sys.stderr.write("unknown file (.NODE) content! - line "+str(c)+"\n")
			exit(1)
		else:
			nodeList2.append([float(cont.group(1)), float(cont.group(2)), float(cont.group(3))])
"""
inputfile.close()



# read the input TOP file
inputfile = open("in/input.top",'r')

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

xxx = """
#write - MESH
outputfile = open ("mesh.vtk", 'w')
	#header
outputfile.write("# vtk DataFile Version 3.1\n"+
				"*description-name of input file*\n"+
				"ASCII\nDATASET UNSTRUCTURED_GRID\n"+
				"POINTS "+str(len(nodeList1)+len(nodeList2))+" FLOAT\n")
	#points
for i in nodeList1:
	outputfile.write(" ".join(i)+"\n")
for i in nodeList2:
	outputfile.write(" ".join(i)+"\n")

	#cells
outputfile.write("\nCELLS "+str(len(elementList1))+" "+str(5*len(elementList1))+"\n")
for i in elementList1:
	outputfile.write("4 "+" ".join(i)+"\n")

	#cell_types
outputfile.write("\nCELL_TYPES "+str(len(elementList1))+"\n")
ss = []
for i in range(len(elementList1)):
	ss.append("10")
outputfile.write(" ".join(ss)+"\n")
outputfile.close()
"""

################################################################################
## input:	.node, .surf
## output:	mesh_surf
################################################################################
# read the input SURF file
inputfile = open("in/input.surf",'r')

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



xxx = """
#write - MESH_SURFACE
outputfile = open ("mesh_surface.vtk", 'w')
	#header
outputfile.write("# vtk DataFile Version 3.1\n"+
				"*description-name of input file*\n"+
				"ASCII\nDATASET UNSTRUCTURED_GRID\n"+
				"POINTS "+str(len(nodeList1))+" FLOAT\n")
	#points - It is not a mistake, list of points AGAIN!
for i in nodeList1:
	outputfile.write(" ".join(i)+"\n")

	#cells
outputfile.write("\nCELLS "+str(len(surfList))+" "+str(4*len(surfList))+"\n")
for i in surfList:
	outputfile.write("3 "+" ".join(i)+"\n")

	#cell_types
outputfile.write("\nCELL_TYPES "+str(len(surfList))+"\n")
ss = []
for i in range(len(surfList)):
	ss.append("5")
outputfile.write(" ".join(ss)+"\n")
outputfile.close()
"""
################################################################################
#prepare description of mesh structure for write

meshStr = prepareMeshStructure(elementList1, 4, "10")
meshSurfStr = prepareMeshStructure(surfList, 3, "5")

writeMesh("out/mesh",0, meshStr, nodeList1, nodeList2) #write mesh
writeMesh("out/meshSurf",0, meshSurfStr, nodeList1,[]) #write meshSurf
exit()

################################################################################
## input:	trajectory.out
## output:	?
################################################################################
# read the input trajectory.out file
inputfile = open("in/trajectory.out",'r')

state = 0
nodeMove = {} #surface_faces

blobID = -1
step = 0

c=0
for line in inputfile.readlines():
	c = c+1

	if state == 1: #DYNAMIC / STATIC info
		state = 2
	elif state == 2: #list of node moves
		cont = re.match(r"(-?\d+.\d+)(e[-\+]\d+)? (-?\d+.\d+)(e[-\+]\d+)? (-?\d+.\d+)(e[-\+]\d+)?\s*",line)
# >>> 3.85e-2
# 0.0385
# >>> 3.85*1e-2
# 0.0385
		if cont == None: #another blob - next step
			state = 0
		else:
			ls = []
			for i in [2,4,6]:
				if cont.group(i) != None:
					ls.append(float(cont.group(i-1))*float("1"+cont.group(i)))
				else:
					ls.append(float(cont.group(i-1)))
			nodeMove[blobID].append(ls)
	if state == 0:	#waiting for blob info
		cont = re.match(r"Blob (\d+), Conformation \d+, step (\d+).*",line)
		if cont != None:
			state = 1
			stepRead = int(cont.group(2))
			if step != stepRead:								
				#write
				#actualize node position = +trajectory
				for i in range(len(nodeList1)):
					for j in range(len(nodeList1[i])):
						nodeList1[i][j] += nodeMove[blobID][i][j]


				for i in range(len(nodeList2)):
					for j in range(len(nodeList2[i])):
						nodeList2[i][j] += nodeMove[blobID][i+len(nodeList1)][j]

				
				writeMesh("out/"+str(blobID)+"mesh",step, meshStr, nodeList1, nodeList2) #write mesh
				writeMesh("out/"+str(blobID)+"meshSurf",step, meshSurfStr, nodeList1,[]) #write meshSurf

				step = stepRead
#				exit()
			blobID = int(cont.group(1))
			nodeMove[blobID] = []
inputfile.close()


















