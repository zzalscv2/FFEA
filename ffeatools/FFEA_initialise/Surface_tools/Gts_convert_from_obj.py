import sys, os
if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + "[INPUT input_fname (.obj)] [OUTPUT output_fname (.gts)]\n")
	
ifile = sys.argv[1]
ofile = sys.argv[2]
verts = 0
edges = 0
faces = 0

sta = open(ifile,'r')
stb = open(ofile,'w')
for i in sta: 
  if i[:2] == "v ":
    verts += 1 
  elif i[:2] == "f ":
    faces += 1  
    edges += 3

# go back to the top,
sta.seek(0)
#   and write the vertices:
stb.write(str(verts) + ' ' + str(edges) + ' ' + str(faces) + ' GtsSurface GtsFace GtsEdge GtsVertex\n')
for i in sta: 
  if i[:2] == "v ":
    stb.write(i[2:])

# go back to the top,
sta.seek(0)
#   and write triangular stuff:
for i in sta:
  if i[:2] == "f ":
    l = i.split()
    v1 = str(int(l[1].split("//")[0]))
    v2 = str(int(l[2].split("//")[0]))
    v3 = str(int(l[3].split("//")[0]))
    stb.write(v1 + " " + v2 + "\n")
    stb.write(v2 + " " + v3 + "\n")
    stb.write(v3 + " " + v1 + "\n")

sta.close()

##### and finally: 
for i in range(faces):
  e1 = 3*i + 1 
  stb.write(str(e1) + " " + str(e1 + 1) + " " + str(e1 + 2) + "\n")

stb.close()
