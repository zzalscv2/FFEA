# -*- coding: utf-8 -*-
# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

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
    """
    The FFEA_turbotrajectory class is an alternative to the FFEA_trajectory
    class. It can only be generated from an FFEA trajectory, and the process
    is one-way. It is not intended for general use, but may be handy for
    applications in which the performance matters more than the semantics.
    In addition to being very fast, it's also very small, as it stores
    no second-order nodes, and stores trajectories as binary blobs, not as
    text files.
    For this reason, the read speed is only limited by your hard drive speed, 
    unlike FFEA_trajectory, which is limited by how fast the individual lines
    from the file can be read in.
    FFEA_turbotrajectory is only really userful in situations when a trajectory
    object has to be read in more than once - it still reads from the .ftj file,
    so ultimately, you're adding a bit of overhead.
    It can also directly precalculate the calls to CGL that are needed for the
    PyMOL plugin to draw triangles on the screen. The major overhead in loading
    trajectories into the viewer is actually looking up the triangle indices
    and calculating the surface normals, and the cgo objects have that
    pre-generated.
    In the future, the FFEA runner may be able to write directly to turbotraj-
    compatible binary blobs.
    """

    def __init__(self):
        self.dimensions = 3 # well ya never know
        return
    
    def load_traj(self, path): # convert trajectory object to turbotraj
        """
        Convert a trajectory into a turbotrajectory.
        Note: this is kinda slow, you're better off using the
        populate_turbotraj_from_ftj function.
        In: self, path to trajectory file.
        Out: writes to the self.turbotraj object.
        """
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
        """
        Get the normal to a plane, specified by three nodes.
        In: node1, node2, node3.
        Out: Plane normal vector as a python list.
        """
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
        In: self, script (an FFEA script object) and display_params, which is
        a dictionary created by the pymol viewer gui.
        Out: populates the self.cgo and self.cgo.blob_index objects.
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
        """
        This writes to the CGO object using a list of elements supplied by the
        user. It outputs a separate pymol object containing only the elements
        that the user specifies. IN order to do this, it has to look up the
        nodes in each element, and then create four triangles for each.
        In: self, element_list (a python list object containing element indices),
        a script object, and a turbotraj object.
        Out: Appends to the cgo.
        """
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
        """            If we had to do this for every node, it would probably be a pain
        Load header info. This works in almost the exact same way as the regular
        FFEA trajectory method.
        """
    
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
        """
        Create a turbotraj object straight from a ,ftj file. As far as I know,
        this is the fastest way to do it, but if you have a better idea, go
        ahead! We load in the file, grab the header, use the first 3 frames
        to get the number of steps (it only starts stepping properly between
        frame 2 and 3) in each frame (why that's not in the header I have no
        idea). 
        From there, we start reading the file sequentially. Every frame of every
        blob has what is essentially header information above it, and let me
        tell ya, it's a lot more human readable than it is machine-readable.
        Regex is used to work out which part of the header we're in, and
        extract the data from the header. We populate the turbotraj object
        only with the first-order nodes.
        In: self, fname...
        Out: populates turbotraj attribute.
        """
        def skip_line(ftj, n):
            for i in range(n):
                ftj.readline()
                
        def match_line(line):
            """
            Match a line according to the highly sophisticated regex query:
            "[0-9]". Uses python's built in regex library. We have to do a bit
            of string fiddling beforehand, because the library takes input and
            gives output in an odd way.
            In: a line (string
            Out: the contents of the line (we use it for blob, conf and step)
            """
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
            """
            Because the .trj file format doesn't have any way of telling you
            what frame you're on, you have to either keep track of it in a
            variable (which gets a bit ungainly because then you also have to
            keep track of which blob you're in) OR just infer it from the steps.
            If you know how many steps there are in each frame (which we do)
            then we can work out which frame we're on.
            """
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
            if line.startswith('Blob 0') and "->" not in line and "Nodes" not in line: # Only lines with blob, conf, step
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
                
                    #try:
                    #    turbotraj[blob][conf][frame][node] = np.fromstring(line, dtype=float, sep=' ')[0:3]
                    node_line = ftj.readline().split()
                    for i in [0,1,2]: # 3 points
                        turbotraj[blob][conf][frame][node][i] = float(node_line[i])
                    # fill the contents of the node with the first 3 numbers in the line, converted from a string. first 3 numbers = ignore second order elements
            if line == "" or line == None: #empty string (no /n) at eof
                break

        self.ftj =  ftj
        self.turbotraj = turbotraj
                    
        print("Trajectory loaded.")
                
                
            
    
class _cgo:
    """
    These are the builtin constants for the PyMOL library's cgo, stolen from
    the source code from that library itself. Reason being, if you're not in
    a pymol extension, you can't actually import PyMOL.
    """
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