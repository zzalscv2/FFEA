#!/usr/bin/env python
import os, sys

if len(sys.argv) != 9:
	sys.exit("Usage: ./traj_to_dx.py [INPUT TRAJECTORY FNAME] [BLOB NUMBER] [NODE NUMBER] [NX] [NY] [NZ] [SCALING] [OUTPUT DX MAP FNAME]")

traj_fname = sys.argv[1]
blob_no = sys.argv[2]
node_no = sys.argv[3]
nx = sys.argv[4]
ny = sys.argv[5]
nz = sys.argv[6]
scaling = sys.argv[7]
dx_fname = sys.argv[8]

temp_fname = ".__temp__" + traj_fname + "_" + dx_fname + "__temp__"

path = os.path.dirname(sys.argv[0])
if path == '':
        path = '.'

os.system(path + "/extract " + traj_fname + " " + blob_no + " " + node_no + " " + temp_fname + "\n")
os.system("python " + path + "/make_dx_map.py " + temp_fname + " " + nx + " " + ny + " " + nz + " " + scaling + " " + dx_fname + "\n")
