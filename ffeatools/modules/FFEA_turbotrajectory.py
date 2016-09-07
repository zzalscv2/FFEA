# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:22:28 2016

@author: rob
"""
import numpy as np
import FFEA_trajectory

class FFEA_turbotrajectory:

    def __init__(self):
        self.dimensions = 3
        return
    
    def load_traj(self, path):
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

    def dump(self, path):
        np.save(path, self.turbotraj)
    
    def pretend_blob(blob):
        print("nothin' here yet")
        return
                