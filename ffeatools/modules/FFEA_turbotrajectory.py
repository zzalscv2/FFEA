# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:22:28 2016

@author: rob
"""
import numpy as np
import FFEA_trajectory
#from pymol import cmd
from pymol.cgo import *
import pymol.cgo as _cgo


class FFEA_turbotrajectory:

    def __init__(self):
        self.dimensions = 3
        return
    
    def load_traj(self, path):
        self.path = path.split(".")[0]
        if path.endswith("out"):
            traj = FFEA_trajectory.FFEA_trajectory(path)
            self.blob = traj.blob # backward compatible
            blobs = traj.blob
            max_blobs = len(traj.blob)
            
            num_confs = []
            for blob in blobs:
                num_confs.append(len(blob))
            
            max_confs = max(num_confs)
            
            num_frames = len(blobs[0][0].frame)
            num_nodes = len(blobs[0][0].frame[0].pos)

            turbotraj = np.empty([max_blobs, max_confs, num_frames, num_nodes, self.dimensions])
            #later , this will automatically reduce the number of dimensions if the blobs or confs are 1
            #but i don't think it has any speed effect
            
            for i in range(max_blobs):
                for j in range(max_confs):
                    for k in range(num_frames):
                        print(str(int(((float(k)/float(num_frames))/(max_blobs*max_confs))*100))+"% converted")
                        for l in range(num_nodes):
                            for m in range(self.dimensions):
                                turbotraj[i][j][k][l][m] = blobs[i][j].frame[k].pos[l][m]
            self.turbotraj = turbotraj
            return
        elif path.endswith("npy"):
            self.turbotraj=np.load(path)
            return
        else:
            raise IOError("File must have .out (traj) or .npy (turbotraj) extension.")

    def dump_traj(self):
        np.save(self.path, self.turbotraj)
    
    def pretend_blob(blob):
        print("nothin' here yet")
        return
        
    def get_normal(self, node0, node1, node2):
        ax = node1[0] - node0[0]
        ay = node1[1] - node0[1]
        az = node1[2] - node0[2]
        bx = node2[0] - node1[0]
        by = node2[1] - node1[1]
        bz = node2[2] - node1[2]

        return [az * by - ay * bz, ax * bz - az * bx, ay * bx - ax * by]

    def create_cgo(self, script):
    """
    This sript creates a cgo object that can be fed into pymol for each frame
    and blob. The cgo object contains the vertices of each triangle in the mesh,
    a normal vector, and some control codes. This method also creates an index
    file listing the name of each cgo object and which frame it belongs to.
    """
      
        def setup(self):
            frames = range(len(self.turbotraj[0][0]))
            surfs = []
            cgo = [] #  this is the cgo list! The _cgo object is the constants.
            cgo_blob_index = []
    
            # cerate a list of surfaces, one for each blob
            for i in range(len(self.turbotraj)):
                surfs.append(script.load_surface(i)) 
            return surfs, frames, cgo, cgo_blob_index
            
        def get_nodes_in_face(turbotraj, face):
            return [turbotraj[blob_num][0][frame][face.n[0]], turbotraj[blob_num][0][frame][face.n[1]], turbotraj[blob_num][0][frame][face.n[2]]]
        
        surfs, frames, self.cgo, self.cgo_blob_index = setup(self)
    
        # for every frame, create a cgo object
        for frame in frames:
            print("Creating frame "+str(frame)+"...")
            sol = [ _cgo.BEGIN, _cgo.TRIANGLES ]
    
            # for each face in each surf, load the nodes into the cgo as triangles
            for blob_num in range(len(surfs)):
                for face in surfs[blob_num].face:
                    nodexyz = get_nodes_in_face(self.turbotraj, face)
                    norm = self.get_normal(nodexyz[0], nodexyz[1], nodexyz[2])
                    sol.extend( [ _cgo.NORMAL, -norm[0], -norm[1], -norm[2], _cgo.VERTEX, nodexyz[0][0]*1000000000, nodexyz[0][1]*1000000000, nodexyz[0][2]*1000000000, _cgo.VERTEX, nodexyz[1][0]*1000000000, nodexyz[1][1]*1000000000, nodexyz[1][2]*1000000000, _cgo.VERTEX, nodexyz[2][0]*1000000000, nodexyz[2][1]*1000000000, nodexyz[2][2]*1000000000 ] )
            sol.append(_cgo.END)
            self.cgo.append(sol)
            # save it for later
            self.cgo_blob_index.append(["blob_"+str(blob_num), frame])

    def dump_cgo(self):
        cgo_array = np.array(self.cgo)
        cgo_blob_index_array = np.array(self.cgo_blob_index)
        np.save(self.path+"_cgo", cgo_array)
        np.save(self.path+"_cgoindex", cgo_blob_index_array)