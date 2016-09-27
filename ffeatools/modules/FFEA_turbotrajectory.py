# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:22:28 2016

@author: rob
"""
import numpy as np
import FFEA_trajectory
#from pymol import cmd
#from pymol.cgo import *
#import pymol.cgo as _cgo
from os import path
import sys
import re as _re

class FFEA_turbotrajectory:

    def __init__(self):
        self.dimensions = 3
        return
    
    def load_traj(self, path): # convert trajectory object to turbotraj
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

    def create_cgo(self, script, display_params):
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

        # optional params
        if display_params['highlight'] != '':
            print("Highlighted nodes "+display_params['highlight']+" detected")
            highlighted = np.fromstring(display_params['highlight'], sep=",", dtype=int) 
            self.write_cgo_for_elements(highlighted, script, self.turbotraj)

    def write_cgo_for_elements(self, element_list, script, turbotraj):
        frames = xrange(len(self.turbotraj[0][0]))
        tops = []
        for i in xrange(len(turbotraj)):
            tops.append(script.load_topology(i)) 
        for frame in frames:
            print("Highlighting nodes in frame "+str(frame))
            sol = [ _cgo.BEGIN, _cgo.TRIANGLES ]
            for blob_num in xrange(len(tops)):
                for element in element_list:
                    nodes = tops[blob_num].element[element].n[:4]
                    face1 = [nodes[0], nodes[1], nodes[2]]
                    face2 = [nodes[1], nodes[2], nodes[3]]
                    face3 = [nodes[0], nodes[2], nodes[3]]
                    face4 = [nodes[0], nodes[1], nodes[3]]
                    faces = [face1, face2, face3, face4]
                    for face in faces:
                        nodexyz = [turbotraj[blob_num][0][frame][face[0]], turbotraj[blob_num][0][frame][face[1]], turbotraj[blob_num][0][frame][face[2]]]
                        norm = self.get_normal(nodexyz[0], nodexyz[1], nodexyz[2])
                        sol.extend( [ _cgo.NORMAL, -norm[0], -norm[1], -norm[2], _cgo.VERTEX, nodexyz[0][0]*1000000000, nodexyz[0][1]*1000000000, nodexyz[0][2]*1000000000, _cgo.VERTEX, nodexyz[1][0]*1000000000, nodexyz[1][1]*1000000000, nodexyz[1][2]*1000000000, _cgo.VERTEX, nodexyz[2][0]*1000000000, nodexyz[2][1]*1000000000, nodexyz[2][2]*1000000000 ] )
            sol.append(_cgo.END)
            self.cgo.append(sol)
            self.cgo_blob_index.append(["highlight", frame])

    def dump_cgo(self):
        cgo_array = np.array(self.cgo)
        cgo_blob_index_array = np.array(self.cgo_blob_index)
        np.save(self.path+"_cgo", cgo_array)
        np.save(self.path+"_cgoindex", cgo_blob_index_array)

    def load_ftj_header(self, fname): # load header data from ftj file and create empty turbotraj
    
        self.path = fname.split(".")[0]
            
        # Get a file object and store it
        try:
            self.ftj = open(fname, "r")

        except(IOError):
            raise IOError("\tFailed to open '" + fname + "' for reading.")

        # Now, read only the information from the top of the file

        # Title
        line = self.ftj.readline().strip()
        if line != "FFEA_trajectory_file":
            raise IOError("\tExpected to read 'FFEA_trajectory_file' but read '" + line + "'. This may not be an FFEA trajectory file.")

        self.ftj.readline()
        self.ftj.readline()

        # num_blobs
        try:
            #self.blob = [[FFEA_trajectory.FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

            line = self.ftj.readline()
            self.num_blobs = int(line.split()[3])

        except(IndexError, ValueError):
            raise IOError("\tExpected to read 'Number of Blobs %d' but read '" + line + "'.")

        # num_conformations
        try:
            line = self.ftj.readline()
            sline = line.split()[3:]
            self.num_conformations = [int(s) for s in sline]

        except(IndexError, ValueError):
            raise IOError("\tExpected to read 'Number of Conformations %d %d ....%d' but read '" + line + "'.")

        # num_nodes
        self.num_nodes = [[0 for i in range(self.num_conformations[j])] for j in range(self.num_blobs)]
        for i in range(self.num_blobs):
            try:
                line = self.ftj.readline()
                sline = line.split()[2:]

                for j in range(self.num_conformations[i]):
                    self.num_nodes[i][j] = int(sline[4 * j + 3])

            except(IndexError, ValueError):
                raise IOError("\tExpected to read 'Blob " + str(i) + ": Conformation 0 Nodes %d Conformation 1 Nodes %d....Conformation " + str(self.num_conformations[i] - 1) + " Nodes %d' but read '" + line + "'.")

        # final whitespace until '*' and save the file pos
        while(self.ftj.readline().strip() != "*"):
            pass

        self.fpos = self.ftj.tell()
        
        # cheeky way of getting the number of frames
        self.ftj.seek(0)
        asterisks = self.ftj.read().count("*")
        self.num_frames = (asterisks-1)/2 # 1 frame for every 2 asterisks, plus an extra one at the end
        self.num_nodes = self.num_nodes[0][0] # not sure why it's trapped in those lists
        self.num_conformations = max(self.num_conformations) # get max number of conformations
        
        # Finally, build the objects
        self.turbotraj = np.empty([self.num_blobs, self.num_conformations, self.num_frames, self.num_nodes, self.dimensions])
        
    def populate_turbotraj_from_ftj(self, fname, surf=None, load_all=1, frame_rate = 1, num_frames_to_read = 1000000, start = 0): # load the contents of a .ftj trajectory into a turbotraj object
    
        def skip_line(ftj, n):
            for i in range(n):
                ftj.readline()
                
        def match_line(line):
            line_contents = line.split(',')
            line_contents_formatted = []
            for thing in line_contents:
                chars = _re.findall("[0-9]", thing)
                num = "".join(chars)
                line_contents_formatted.append(int(num))
            try:
                return line_contents_formatted[0], line_contents_formatted[1], line_contents_formatted[2] # blob, conf, step
            except IndexError:
                raise IndexError("Expected line contents to contain 3 items, got "+str(line_contents_formatted)+" of length "+str(len(line_contents_formatted)))
            
        def convert_step_to_frame(steps_per_frame, step):
            if step > 1:
                return ((step-1)/steps_per_frame)+1
            else:
                return step

        print("Loading FFEA trajectory file...")

        # Test file exists
        if not path.exists(fname):
            raise IOError("No trajectory found at that location")

        # Header first, for sure
        self.load_ftj_header(fname)

        try:
            ftj = open(fname, "r")
        except(IOError):
            raise IOError("\tFailed to open '" + fname + "' for reading.")
            
        
        # get step - find the 2nd and 3rd frames, get the step difference
        seek = True
        steps = []
        while seek:
            line = ftj.readline()
            if line.startswith('Blob') and "->" not in line and "Nodes" not in line: # Only lines with blob, conf, step
                blob, conf, step = match_line(line) # regex the line to get the values
                steps.append(step)
            if len(steps) >= 3:
                break
            
        self.step = steps[2] - steps[1]
        steps_per_frame = self.step
        
        print("Steps calculated successfully...")
        print("Loading trajectory...")
        
        turbotraj = self.turbotraj #faster to access a local object
            
        ftj.seek(0)
        
        seek = True
        while seek:
            line = ftj.readline()
            if line.startswith('Blob') and "->" not in line and "Nodes" not in line: #Only lines with blob, conf, step
                blob, conf, step = match_line(line) # regex the line to get the values
                frame = convert_step_to_frame(steps_per_frame, step)
                ftj.readline() # skip the word 'DYNAMIC'
                nodes_range = xrange(self.num_nodes)
                for node in nodes_range: #  each node occupies one line
                    #self.turbotraj[blob][conf][frame][node] = np.fromstring(self.ftj.readline(), dtype=float, sep=' ')[0:3]
                    node_line = ftj.readline().split()
                    for i in [0,1,2]: # 3 points
                        turbotraj[blob][conf][frame][node][i] = float(node_line[i])
                    # fill the contents of the node with the first 3 numbers in the line, converted from a string. first 3 numbers = ignore second order elements
            if line == "": #empty string (no /n) at eof
                break

        self.ftj =  ftj
        self.turbotraj = turbotraj
                    
        print("Trajectory loaded.")
                
                
            
    
class _cgo:
    POINTS             = 0.0
    LINES              = 1.0
    LINE_LOOP          = 2.0
    LINE_STRIP         = 3.0
    TRIANGLES          = 4.0
    TRIANGLE_STRIP     = 5.0
    TRIANGLE_FAN       = 6.0
    STOP               =  0.0
    NULL               =  1.0
    BEGIN              =  2.0
    END                =  3.0
    VERTEX             =  4.0
    NORMAL             =  5.0
    COLOR              =  6.0
    SPHERE             =  7.0
    TRIANGLE           =  8.0
    CYLINDER           =  9.0
    LINEWIDTH          = 10.0
    WIDTHSCALE         = 11.0
    ENABLE             = 12.0
    DISABLE            = 13.0
    SAUSAGE            = 14.0
    CUSTOM_CYLINDER    = 15.0
    DOTWIDTH           = 16.0
    ALPHA_TRIANGLE     = 17.0
    ELLIPSOID          = 18.0
    FONT               = 19.0
    FONT_SCALE         = 20.0
    FONT_VERTEX        = 21.0
    FONT_AXES          = 22.0
    CHAR               = 23.0
    ALPHA              = 25.0
    QUADRIC            = 26.0 # NOTE: Only works with ellipsoids and disks
    CONE               = 27.0 
    LIGHTING           = float(0x0B50)