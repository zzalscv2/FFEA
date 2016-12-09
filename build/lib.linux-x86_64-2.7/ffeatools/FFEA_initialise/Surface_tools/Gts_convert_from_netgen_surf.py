import sys, os
if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + "[INPUT input_fname (.surf)] [OUTPUT output_fname (.gts)]\n")
	
ifile = sys.argv[1]
ofile = sys.argv[2]
verts = 0
edges = 0
faces = 0

sta = open(ifile,'r')
stb = open(ofile,'w')
sta.readline()
verts = int(sta.readline())
for i in range(verts):
	sta.readline()
faces = int(sta.readline())
edges = faces * 3

# go back to the top,
sta.seek(0)
#   and write the vertices:
stb.write(str(verts) + ' ' + str(edges) + ' ' + str(faces) + ' GtsSurface GtsFace GtsEdge GtsVertex\n')
sta.readline()
sta.readline()
for i in range(verts): 
	stb.write(sta.readline())

#   and write triangular stuff:
sta.readline()
for i in range(faces):
	l = sta.readline().split()
	v1 = l[0].strip()
	v2 = l[1].strip()
	v3 = l[2].strip()
	stb.write(v1 + " " + v2 + "\n")
	stb.write(v2 + " " + v3 + "\n")
	stb.write(v3 + " " + v1 + "\n")

sta.close()

##### and finally: 
for i in range(faces):
  e1 = 3*i + 1 
  stb.write(str(e1) + " " + str(e1 + 1) + " " + str(e1 + 2) + "\n")

stb.close()
