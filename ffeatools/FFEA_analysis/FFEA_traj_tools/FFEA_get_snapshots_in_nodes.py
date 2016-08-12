import argparse, os 
from FFEA_script import FFEA_script
from FFEA_trajectory import FFEA_trajectory
from FFEA_node import FFEA_node
from math import ceil

def parser():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--ffea', help="ffea script")
  parser.add_argument('-t', '--traj', help="input file", type=argparse.FileType('r'))
  parser.add_argument('-o', '--oName', help="output basename")
  parser.add_argument('-s', '--snapshots', help="list of snapshots required, e. g., 5,12,44-48")
  return parser


def checkArgs(a):
  if (a.ffea == None):
     print "A valid FFEA script is needed"
     return False
  if (os.path.isfile(a.ffea)):
     a.ffea = FFEA_script(a.ffea)
  else: 
     print a.ffea, " did not work as a valid FFEA script"
  if (a.traj == None): 
     a.traj = a.ffea.params.trajectory_out_fname
     print "we'll use ", a.traj, " as the trajectory to be read"
  if (a.oName == None):
     a.oName, crap = os.path.splitext(a.traj)
     print "we'll use ", a.oName, " as a base name for writing the nodes"
  if (a.snapshots == None):
     print "a list of snapshots needs to be specified"
     return False
  else:
     list = []
     try:
       for i in a.snapshots.split(","):
         if i.count("-"):
           for j in range(int(i.split("-")[0]), int(i.split("-")[1]) +1):
             list.append(j)
         else:
           list.append(int(i))
     except:
       print "invalid list of snapshots!", a.snapshots
       return False
     list.sort()
     a.snapshots = list

  return True



p = parser()
args = p.parse_args()
if (checkArgs(args) == False):
  p.print_usage()
  quit()

Traj = FFEA_trajectory(args.traj, load_all=0)
Traj.skip_frame()

cnt = 0
for i in args.snapshots:
  print i
  while cnt != i:
    Traj.skip_frame()
    cnt += 1 

  print cnt 
  Traj.load_frame()
  cnt += 1 


NIN = []
NSN = []
for i in range(args.ffea.params.num_blobs):
  tmpI = []
  tmpS = []
  for j in range(args.ffea.params.num_conformations[i]):
    nodefile = args.ffea.blob[i].conformation[j].nodes
    Node = FFEA_node(nodefile)
    tmpS.append( Node.num_surface_nodes ) 
    tmpI.append( Node.num_interior_nodes ) 
  NIN.append(tmpI)
  NSN.append(tmpS)



for ib in range(len(Traj.blob)):
  for ic in range(len(Traj.blob[ib])):
    for iF in range(len(Traj.blob[ib][ic].frame)):
      Traj.blob[ib][ic].frame[iF].num_interior_nodes = NIN[ib][ic]
      Traj.blob[ib][ic].frame[iF].num_surface_nodes = NSN[ib][ic]
     
    
      if Traj.blob[ib][ic].motion_state == "STATIC":
        print "STATIC"
      else:
        Traj.blob[ib][ic].frame[iF].pos /= args.ffea.blob[ib].scale
        oFile = args.oName + "_blob" + str(ib) + "-conf" + str(ic) + "-step" + str(args.snapshots[iF]) + ".nodes" 
        Traj.blob[ib][ic].frame[iF].write_to_file(oFile)





