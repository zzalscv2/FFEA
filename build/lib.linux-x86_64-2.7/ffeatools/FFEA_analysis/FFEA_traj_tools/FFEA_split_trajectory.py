import argparse, os 
from math import ceil

def parser():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--ifile', help="input file", type=argparse.FileType('r'))
  parser.add_argument('-o', '--ofile', help="output file", type=argparse.FileType('w'))
  parser.add_argument('-b', '--begin', help="first snapshot included", type=int)
  parser.add_argument('-e', '--end', help="last snapshot, included", type=int)
  parser.add_argument('-p', '--parts', help="[OPTIONAL] number of parts", type=int)
  return parser

def checkArgs(a):
  ## ifile and ofile:
  if (a.ifile == None): 
     print "invalid input file: ", a.ifile
     return False
  if (a.ofile == None):
     print "invalid output file: ", a.ofile
     return False 
  if (a.begin == None):
     print "starting at snapshot 0"
     a.begin = 0
  elif (a.begin < 0):
     print "invalid first snapshot!", a.begin
     return False
  if (a.end == None) or (a.end < 0):
     print "invalid last snapshot!", a.end
     return False
  if (a.begin > a.end): 
     print "first snapshot should be before than the last one"
     return False
  if (a.parts == None):
     a.parts = 1
  elif (a.parts < 0): 
     print "the number of parts needs to be positive!"
     return False
  return True

def forget(ifile,n):
  for i in range(n):
    ifile.readline()
  return

def transfer(ifile, ofile, n):
  for i in range(n):
    ofile.write(ifile.readline())
  return

def getAndWriteLine(ifile,ofile):
  l = ifile.readline()
  ofile.write(l)
  return l


p = parser()
args = p.parse_args()
if (checkArgs(args) == False):
  p.print_usage()
  quit()

## skip header:
transfer(args.ifile,args.ofile,3)

# get the number of blobs:

n_b = int(getAndWriteLine(args.ifile,args.ofile).split()[-1])
print "number of blobs", n_b

# get the number of conformations:
list = getAndWriteLine(args.ifile,args.ofile).split()[3:]
if len(list) != n_b:
  print "ABORTING: the number of blobs does not match the conformations list"
  quit()
n_c = []
for i in list:
  n_c.append(int(i))
print "conformations: ", n_c

# get the number of nodes per blob/conformation:
D = {} ## dictionary with the number of nodes. 
       # It has the form (0,2) = 454, i. e., blob 0, conformation 2, has 454 nodes.
for i in range(n_b):
  for j in range(n_c[i]):
    l = getAndWriteLine(args.ifile,args.ofile).split()
    D[(int(l[1][:-1]), int(l[3]))] = int(l[-1])

# finish the header:
transfer(args.ifile, args.ofile, 2)

# setup the parts:
args.ofile.flush()
ofileBaseName = args.ofile.name
H = []
sta = open(ofileBaseName, 'r')
H = sta.readlines()
sta.close()

# read the first conformation:
l = args.ifile.readline()
step = int(l.split()[-1])
# and forget everything until you reach the first snapshot
while (step < args.begin):
  for i in range(n_b):
    c = int(l.split()[3][:-1])
    b = int(l.split()[1][:-1])
    forget(args.ifile, D[b,c]+1)
    l = args.ifile.readline()
  forget(args.ifile, 2 + n_b)
  l = args.ifile.readline()
  step = int(l.split()[-1])

# 
print "Will start writing at step: ", step

chunk = int(ceil((args.end - args.begin)/args.parts))
part = 1
upTo = args.begin + chunk
while (step < args.end):
  if (step > upTo): 
    args.ofile.close()
    if (part == 1): 
      os.rename(ofileBaseName, ofileBaseName + ".part1")
    upTo += chunk
    part += 1
    
    args.ofile = open(ofileBaseName + ".part" + str(part),'w')
    for zz in H:
      args.ofile.write(zz)
 


  args.ofile.write(l)
  for i in range(n_b):
    c = int(l.split()[3][:-1])
    b = int(l.split()[1][:-1])
    transfer(args.ifile, args.ofile, D[b,c]+1)
    l = getAndWriteLine(args.ifile,args.ofile)
  transfer(args.ifile, args.ofile, 2 + n_b)
  
  l = args.ifile.readline()
  try:
    step = int(l.split()[-1])
  except IndexError:
    break

args.ifile.close()
args.ofile.close()




